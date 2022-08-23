#ifndef RAJA_FDConversion_HPP
#define RAJA_FDConversion_HPP

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
        if(sameLayout) { return; }

        T * _v2 = new T[v.size()];
        ViewType v2(_v2, l);
        
        auto bounds = layout_to_bounds<NumDims>(l);
        auto lambda = conversion_lambda<T, LayoutType>(v, v2);

        kernel<kpol>(bounds, lambda);

        T * _v = v.get_data();   
        memcpy(_v, _v2, sizeof(T) *v.size());
        v.set_layout(l);
        delete[] _v2;
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

} //namespace RAJA
#endif