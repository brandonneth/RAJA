#ifndef RAJA_LCAutoLayoutUtils_HPP
#define RAJA_LCAutoLayoutUtils_HPP
#include <ctime>
namespace RAJA {

template <camp::idx_t CurrDim, camp::idx_t NumDims>
struct nest_helper;

template <camp::idx_t NumDims>
struct nest_helper<0,NumDims> {
  using pol = statement::Lambda<0>;
};

template <camp::idx_t CurrDim, camp::idx_t NumDims>
struct nest_helper {
  using subpol = typename nest_helper<CurrDim-1, NumDims>::pol;  
  
  using pol = statement::For<NumDims-CurrDim, omp_parallel_for_exec, subpol>;
};

template <camp::idx_t NumDims>
struct nest {
  using _pol = typename nest_helper<NumDims, NumDims>::pol;
  using pol = KernelPolicy<_pol>;
};

template <camp::idx_t DimsLeft>
auto layout_to_bounds_helper(std::array<camp::idx_t, DimsLeft> sizes) {

  if constexpr (DimsLeft == 0) {
    return make_tuple();
  } else {
    auto thisdim = make_tuple(RangeSegment(0,sizes[0]));

    std::array<camp::idx_t, DimsLeft-1> remaining_sizes; 
    for(int i = 0; i < DimsLeft-1; i++) {
      remaining_sizes[i] = sizes[i+1];
    } 

    auto other_bounds = layout_to_bounds_helper<DimsLeft-1>(remaining_sizes);
    return tuple_cat(thisdim, other_bounds);
  }
}

template <camp::idx_t NumDims>
auto layout_to_bounds(Layout<NumDims> l) {

  std::array<camp::idx_t, NumDims> sizes;
  for(int i = 0; i < NumDims; i++) {
    sizes[i] = l.sizes[i];
  }

  return layout_to_bounds_helper<NumDims>(sizes); 
  
}

template <typename Args, camp::idx_t...Is, typename ViewType>
void cl_helper(Args&& args_tuple, camp::idx_seq<Is...>, ViewType & v, ViewType &v2) {
  auto val = v(camp::get<Is>(args_tuple)...);
  v2(camp::get<Is>(args_tuple)...) = val; 
}

template <typename T, typename LayoutType>
auto conversion_lambda(View<T, LayoutType> & v, View<T, LayoutType> &v2) {
  return [&](auto...is) {
    cl_helper(camp::forward_as_tuple(is...), idx_seq_from_to<0,sizeof...(is)>(), v, v2);
  };
}


template <typename ViewType, typename LayoutType>
auto permute_view(ViewType & v, LayoutType l) {

  using T = typename ViewType::value_type;
  constexpr camp::idx_t NumDims = l.n_dims;
  using kpol = typename nest<NumDims>::pol;

  auto converter = [=, &v](auto) {

    int sameLayout = 1;
    
    for(auto i = 0; i < NumDims; i++) {
      if(v.get_layout().strides[i] != l.strides[i]) {
        sameLayout = 0;
      }
    }
    if(sameLayout) {
      //std::cerr << "no layout change\n"; 
      return;
    }
    T * _v2 = new T[v.size()];

    ViewType v2(_v2, l);
    auto bounds = layout_to_bounds<NumDims>(l);

    auto lambda = conversion_lambda<T, LayoutType>(v, v2);
    
    
    kernel<kpol>(bounds, lambda);
    
    T * _v = v.get_data();   
    memcpy(_v, _v2, sizeof(T) *v.size());
    v.set_layout(l);
    free(_v2);
  };

  auto knl = make_forall<loop_exec>(RangeSegment(0,1), converter);
  return knl;


}


template <typename... ViewTypes, typename... LayoutTypes, camp::idx_t...Is>
auto permute_views(camp::tuple<ViewTypes...> views, camp::tuple<LayoutTypes...> layouts, camp::idx_seq<Is...>) {
  return make_tuple(permute_view(camp::get<Is>(views), camp::get<Is>(layouts))...);
}

template <typename... ViewTypes, typename... LayoutTypes>
auto permute_views(camp::tuple<ViewTypes...> views, camp::tuple<LayoutTypes...> layouts) {
  auto seq = idx_seq_for(views);
  return permute_views(views, layouts, seq);
}

template <typename ViewTuple, typename KernelType1, typename KernelType2>
auto manipulate_layouts(ViewTuple views, KernelType1 knl1, KernelType2 knl2) {

  auto endlayouts = modeling_magic(views, knl1, knl2);

  auto l1 = get<0>(endlayouts);
  auto conversion_knls = permute_views(views, endlayouts);

  auto all_knls = tuple_cat(make_tuple(knl1), conversion_knls, make_tuple(knl2));
  
  auto computation = chain(all_knls);

  return computation; 
}


template <typename T>
auto layouts_and_names(View<T, Layout<2>> view) {

  std::array<camp::idx_t, 2> sizes;

  for(int i = 0; i < 2; i++) {
    sizes[i] = view.get_layout().sizes[i];
  }


  auto layout1 = make_permuted_layout(sizes, {{0,1}});
  auto layout2 = make_permuted_layout(sizes, {{1,0}});

  auto name1 = "0_1_to_0_1";
  auto name2 = "0_1_to_1_0";
  
  return make_tuple(make_tuple(layout1, name1), make_tuple(layout2, name2));

  

}


template <typename ViewType, typename KT1, typename KT2, typename LNTuple>
auto create_variant(ViewType & view, KT1 knl1, KT2 knl2, LNTuple layout_and_name) {

  auto conv_knl = permute_view(view, get<0>(layout_and_name));

  auto knls = make_tuple(knl1, conv_knl, knl2);

  auto computation = chain(knls);

  return make_tuple(computation, get<1>(layout_and_name));

}
template <typename ViewType, typename KT1, typename KT2, typename LNTuple, camp::idx_t...Is>
auto create_variants(ViewType & view, KT1 knl1, KT2 knl2, LNTuple layouts_and_names, camp::idx_seq<Is...>) {

  return make_tuple(create_variant(view, knl1, knl2, get<Is>(layouts_and_names))...);
}


template <typename ViewType, typename KernelType1, typename KernelType2>
auto all_variants_one_view(ViewType & view, KernelType1 knl1, KernelType2 knl2) {

  auto end_layouts_and_names = layouts_and_names(view);
  
  auto variants = create_variants(view, knl1, knl2, end_layouts_and_names, idx_seq_for(end_layouts_and_names));

  return variants;

}

template <typename ViewTuple, typename KernelType1, typename KernelType2>
auto all_variants(ViewTuple views, KernelType1 knl1, KernelType2 knl2) {

  return all_variants_one_view(get<0>(views), knl1, knl2);
  

}


template <typename...Contents>
auto named_kernels_executor_helper(tuple<Contents...> named_kernels) {
  if constexpr (sizeof...(Contents) == 1) {
    auto knl = get<0>(get<0>(named_kernels));
    auto name = get<1>(get<0>(named_kernels));
    auto start = std::clock();
    knl();
    auto stop = std::clock();
    auto t = stop - start;

    std::cout << name << ", " << t << "\n";
    return;
  } else {
    auto first_knl = make_tuple(get<0>(named_kernels));
    auto rest_knls = tuple_slice<1, sizeof...(Contents)>(named_kernels);
     named_kernels_executor_helper(first_knl);
   named_kernels_executor_helper(rest_knls);
  }
}

template <typename...Contents>
auto execute_named_kernels(tuple<Contents...> named_kernels) {
  named_kernels_executor_helper(named_kernels);
}
}//namespace RAJA


#endif

