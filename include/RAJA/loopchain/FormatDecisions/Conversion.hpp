#ifndef RAJA_FDConversion_HPP
#define RAJA_FDConversion_HPP
#include <cmath>
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
  
  using pol = statement::For<NumDims-CurrDim, loop_exec, subpol>;
};

template <camp::idx_t NumDims>
struct nest {
  using _pol = typename nest_helper<NumDims, NumDims>::pol;
  using pol = KernelPolicy<_pol>;
};

template <camp::idx_t DimsLeft>
auto layout_to_bounds_helper(std::array<camp::idx_t, DimsLeft> sizes) {

  if constexpr (DimsLeft == 0) {
    camp::sink(sizes); // quiets use-def warning
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

template <typename Args, camp::idx_t...Is, typename ViewType>
void c_helper(Args&& args_tuple, camp::idx_seq<Is...>, ViewType & v, ViewType &  v2) {
  v(camp::get<Is>(args_tuple)...) = v2(camp::get<Is>(args_tuple)...);

}

template <typename ViewType>
auto c_lambda(ViewType & v, ViewType & v2) {
  return [&](auto...is) {
    c_helper(camp::forward_as_tuple(is...), idx_seq_from_to<0,sizeof...(is)>(), v, v2);
  };
}
template <camp::idx_t I>
camp::idx_t estimate_conversion_cost(std::array<camp::idx_t, I> in_perm, std::array<camp::idx_t, I> out_perm) {
  static std::map<std::string,idx_t> cache;
  std::stringstream st;
  for(size_t i = 0; i < in_perm.size(); i++) {
    st << in_perm[i] << " " << out_perm[i];
  }
  std::string key = st.str();
  if(cache[key] != 0) {
    return cache[key];
  }
  if (in_perm == out_perm) {
    cache[key] = 1;
    return cache[key];
  }

  using V =  View<double, Layout<I>>;
  constexpr camp::idx_t problem_size = 20000;
 double root = std::pow(problem_size, 1.0/I);
 camp::idx_t n = (camp::idx_t) root;

  camp::idx_t s = 1;
  std::array<camp::idx_t, I> sizes;
  for(camp::idx_t i = 0; i < I; i++) {
    sizes[i] = n;
    s *= n;
  }
  Layout<I> in_layout = make_permuted_layout(sizes, in_perm);
  Layout<I> out_layout = make_permuted_layout(sizes, out_perm);
  
  using kpol = typename nest<I>::pol;

  auto start =std::chrono::high_resolution_clock::now(); 
  V in(new double[s], in_layout);

  V out(new double[s], out_layout);

  auto l = c_lambda(in, out);
  auto seg = RangeSegment(0,n);

  auto segs = tuple_repeat<I>(seg);

  auto conversion = make_kernel<kpol>(segs, l);

  conversion(); 
  using T = double;
T * _v = in.get_data();   
  T * _v2 = out.get_data();
  memcpy(_v, _v2, sizeof(double) *in.size());
    free(_v2);
    free(_v);
  auto stop =std::chrono::high_resolution_clock::now();

  auto t =std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start).count();
  cache[key] = t;
  //std::cerr << "convcost " << key << ": " << t << "\n";
  return cache[key];

} 

} //namespace RAJA
#endif
