// Contains code that implements the fuse transformation for kernels
// Because this transformation does not analyze dependences and just applie
// the transformation, this file does not include any ISL analysis code
#ifndef RAJA_LCFUSE_HPP
#define RAJA_LCFUSE_HPP

#include "RAJA/config.hpp"
#include "RAJA/loopchain/KernelWrapper.hpp"
#include "RAJA/loopchain/Chain.hpp"
#include "RAJA/loopchain/transformations/Common.hpp"
#include "RAJA/loopchain/transformations/BoundaryKernels.hpp"

namespace RAJA
{

//Transformation Declarations

template <typename... KernelTypes>
auto fuse(camp::tuple<KernelTypes...> knlTuple);

template <typename...KernelTypes>
auto fuse(KernelTypes... knls);


//Transformation Definitions

// returns a lambda that executes the lambdas in lambdas one at a time, in order.
template <typename...LambdaTypes, camp::idx_t...Is>
auto fused_lambda(camp::tuple<LambdaTypes...> lambdas, camp::idx_seq<Is...>) {
  return [=](auto...is) {
    camp::sink((camp::get<Is>(lambdas)(is...), 0)...);
  };
}

// returns a kernel object that executes the kernels in knlTuple in a fused schedule for the 
// parts of their iteration spaces that they all share. For each iteration point in the 
// shared iteration space, the produced kernel will execute the first kernel in the tuple for
// that iteration point, then the second kernel, etc. After each kernel has been executed for 
// that iteration point, the produced kernel moves to the second iteration point and repeats
template <typename...KernelTypes, camp::idx_t...Is>
auto fused_knl(camp::tuple<KernelTypes...> knlTuple, camp::idx_seq<Is...> seq) {

  auto segmentTuples = make_tuple(camp::get<Is>(knlTuple).segments...);
  auto sharedIterSpaceTuple = intersect_segment_tuples(segmentTuples);

  auto lambdaTuple = make_tuple(camp::get<0>(camp::get<Is>(knlTuple).bodies)...);
  auto fusedLambda = fused_lambda(lambdaTuple, seq);
  
  auto sampleKnl = camp::get<0>(knlTuple);
  using KPol = typename decltype(sampleKnl)::KPol;

  return make_kernel<KPol>(sharedIterSpaceTuple, fusedLambda);
}


template <typename... KernelTypes, camp::idx_t...Is>
auto fuse(camp::tuple<KernelTypes...> knlTuple, camp::idx_seq<Is...> seq) {
  auto lowBoundKnls = low_boundary_knls(knlTuple, seq);
  auto highBoundKnls = high_boundary_knls(knlTuple, seq);
  auto fusedKnl = fused_knl(knlTuple, seq);

  auto allKnls = tuple_cat(lowBoundKnls, make_tuple(fusedKnl), highBoundKnls);
  
  return grouped_kernels(allKnls);
}


template <typename... KernelTypes>
auto fuse(camp::tuple<KernelTypes...> knlTuple) {
  return fuse(knlTuple, idx_seq_for(knlTuple));
}

template <typename...KernelTypes>
auto fuse(KernelTypes... knls) {
  return fuse(make_tuple(knls...));
}

} //namespace RAJA

#endif
