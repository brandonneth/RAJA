#ifndef _RAJA_SPARSIFY_HPP
#define _RAJA_SPARSIFY_HPP

namespace RAJA {

std::vector<camp::idx_t> invert_permutation_unsafe(std::vector<camp::idx_t> p) {
  std::vector<camp::idx_t> q = std::vector<camp::idx_t>(p.size());
  for(size_t i = 0; i < p.size(); i++) {
    q[p[i]] = i;
  }
  return q;
}

std::vector<camp::idx_t> apply_permutation(std::vector<camp::idx_t> p, std::vector<camp::idx_t> v) {
  std::vector<camp::idx_t> pp = std::vector<camp::idx_t>(v.size());

  for(size_t i = 0; i < p.size(); i++) {
    pp[p[i]] = v[i];
  }

  return pp;
}

template <typename ElmType, std::size_t ViewDims, std::size_t CurrIdx, std::size_t MaxIdx>
auto sparsification_follow_dimensions( SparseView<ElmType,ViewDims> & A, std::vector<camp::idx_t> argumentOrderInverse, LeadSegmentImpl * lead) {
  if constexpr (CurrIdx == MaxIdx) {
    camp::sink(A,argumentOrderInverse,lead);
    return make_tuple();
  } else {
    
    auto currDimension = make_follow_segment(A.impl->dim(argumentOrderInverse[CurrIdx]), lead);
    auto remainingDims = sparsification_follow_dimensions<ElmType, ViewDims, CurrIdx+1, MaxIdx>(A, argumentOrderInverse, lead);
    return tuple_cat(camp::tuple(currDimension), remainingDims);
  }
}

template <typename KernelPolicy, std::size_t ViewDims, typename ElmType,typename... DenseTupleTypes>
auto sparsify_iteration_space_equal_dims(camp::tuple<DenseTupleTypes...> denseSegs, SparseView<ElmType,ViewDims> & A, std::vector<camp::idx_t> argumentOrder) {
  camp::sink(denseSegs);
  if constexpr (ViewDims != sizeof...(DenseTupleTypes)) {
    std::cerr << "Warning in sparsify_iteration_space_equal_dims: Iteration space and View dimensionality mismatch.\n";
  }

  auto nestingOrder = policy_to_nesting_order<KernelPolicy>();
  auto nestingOrderInverse = invert_permutation_unsafe(nestingOrder);

  auto q = invert_permutation_unsafe(argumentOrder);
 
  auto applied = apply_permutation(nestingOrderInverse, q);
  A.impl->set_permutation(applied); 
  auto leadImpl = new LeadSegmentImpl(A.impl->dim(q[0]));

  auto followDimensions = sparsification_follow_dimensions<double,ViewDims,1,ViewDims>(A, q, leadImpl);

  auto leadDimension = new SparseSegment(leadImpl);
  auto sparseDimensions = tuple_cat(camp::tuple(*leadDimension), followDimensions);

  return sparseDimensions;
  
}

template <typename KernelPolicy, std::size_t ViewDims, typename ElmType,typename... DenseTupleTypes>
auto sparsify_iteration_space_smaller_data(camp::tuple<DenseTupleTypes...> denseSegs, SparseView<ElmType,ViewDims> & A, std::vector<camp::idx_t> argumentOrder) {
  camp::sink(A, argumentOrder);
  if constexpr (ViewDims >= sizeof...(DenseTupleTypes)) {
    std::cerr << "Warning in sparsify_iteration_space_smaller_data:  View dimensionality greater than or equal to iteration space dimensionality.\n";
  }

  auto nestingOrder = policy_to_nesting_order<KernelPolicy>();
  auto nestingOrderInverse = invert_permutation_unsafe(nestingOrder);

  camp::sink(nestingOrderInverse);
  return denseSegs;
}

template <typename KernelPolicy, std::size_t ViewDims, typename ElmType,typename... DenseTupleTypes>
auto sparsify_iteration_space_larger_data(camp::tuple<DenseTupleTypes...> denseSegs, SparseView<ElmType,ViewDims> & A, std::vector<camp::idx_t> argumentOrder) {
  camp::sink(A, argumentOrder);
  if constexpr (ViewDims <= sizeof...(DenseTupleTypes)) {
    std::cerr << "Warning in sparsify_iteration_space_larger_data:  View dimensionality less than or equal to iteration space dimensionality.\n";
  }

  auto nestingOrder = policy_to_nesting_order<KernelPolicy>();
  auto nestingOrderInverse = invert_permutation_unsafe(nestingOrder);

  camp::sink(nestingOrderInverse);
  return denseSegs;
}


template <typename KernelPolicy, std::size_t ViewDims, typename ElmType,typename... DenseTupleTypes>
auto sparsify_iteration_space(camp::tuple<DenseTupleTypes...> denseSegs, SparseView<ElmType,ViewDims> & A, std::vector<camp::idx_t> argumentOrder) {
  if constexpr (ViewDims == sizeof...(DenseTupleTypes)) {
    return sparsify_iteration_space_equal_dims<KernelPolicy>(denseSegs, A, argumentOrder);
  } else if constexpr (ViewDims < sizeof...(DenseTupleTypes)) {
    return sparsify_iteration_space_smaller_data<KernelPolicy>(denseSegs, A, argumentOrder);
  } else {
    return sparsify_iteration_space_larger_data<KernelPolicy>(denseSegs, A, argumentOrder);

    
  
  }
}




}


#endif
