// Contains code that implements the fuse transformation for kernels
// Because this transformation does not analyze dependences and just applie
// the transformation, this file does not include any ISL analysis code
#ifndef RAJA_LCFUSE_HPP
#define RAJA_LCFUSE_HPP

#include "RAJA/config.hpp"
#include "RAJA/loopchain/KernelWrapper.hpp"
#include "RAJA/loopchain/Chain.hpp"
#include "RAJA/loopchain/transformations/Common.hpp"
namespace RAJA
{

//Transformation Declarations

template <typename... KernelTypes>
auto fuse(camp::tuple<KernelTypes...> knlTuple);

template <typename...KernelTypes>
auto fuse(KernelTypes... knls);


//Transformation Definitions


template <typename SegmentTuple, camp::idx_t I>
RAJA_INLINE
auto low_boundary_start_dims_helper(SegmentTuple knlSegments, SegmentTuple sharedSegments) {
  auto knlSegment = camp::get<I>(knlSegments);
  auto sharedSegment = camp::get<I>(sharedSegments);

  using RangeType = decltype(knlSegment);
  return RangeType(*knlSegment.begin(), *sharedSegment.end());
}

template <typename SegmentTuple, camp::idx_t... Is>
RAJA_INLINE
auto low_boundary_start_dims(SegmentTuple knlSegments, SegmentTuple sharedSegments, camp::idx_seq<Is...>) {
  return make_tuple((low_boundary_start_dims_helper<Is>(knlSegments, sharedSegments))...);
}

// Creates one of the hyper-rectangles in the lower boundary 
// iteration space decomposition for a single kernel
<typename...SegmentTypes, camp::idx_t I, camp::idx_t...Is>
auto low_boundary_decomp_piece(camp::tuple<SegmentTypes...> knlSegments, 
                               camp::tuple<SegmentTypes...> sharedSegments,
                               camp::idx_seq<Is...> seq) {
  //For the dimensions up to the piece number, we start at the shared start and go to the original end
  auto startSeq = idx_seq_from_to<0,I>();
  auto startDims = low_boundary_start_dims(knlSegments, sharedSegments, firstSeq);

  //For the dimension at the piece number, start at the original start and go to the start of the shared
  auto ithKnlSegment = camp::get<I>(knlSegments);
  auto ithSharedSegment = camp::get<I>(sharedSegments);
  auto ithDim = make_tuple(RangeSegment(*ithKnlSegment.begin(), *ithSharedSegment.begin()));

  //For the rest of the dimensions its original.start to original.end
  auto endSeq = idx_seq_from_to<I+1, sizeof...(Is)>();
  auto endDims = tuple_slice<I+1, sizeof...(Is)>(knlSegments);

  return tuple_cat(startDims, ithDim, endDims);
} 


// Decomposes the lower boundary iteration space of a single kernel into a 
// tuple of hyper-rectangles. It does so by using the lower faces of the sharedSegment
// as partition hyperplanes
<typename...SegmentTypes, camp::idx_t...Is>
auto low_boundary_decomp(camp::tuple<SegmentTypes...> knlSegments, 
                         camp::tuple<SegmentTypes...> sharedSegments,
                         camp::idx_seq<Is...> seq) {
  return make_tuple(low_boundary_decomp_piece<Is>(knlSegments, sharedSegments, seq));
}


// Returns the kernels that execute the lower boundary iterations for a single
// loop in a loop chain. 
template <typename KPol, typename SegmentTuple, typename Bodies..., typename...SegmentTypes, camp::idx_t...Is>
auto low_boundary_knls_for_knl(KernelWrapper<KPol,SegmentTuple,Bodies...> knl,
                               camp::tuple<SegmenTypes...> sharedSegmentTuple, 
                               camp::idx_seq<Is...> seq) {
  auto lowBoundaryIterSpaceDecomp = low_boundary_decomp(knl.segments, sharedSegmentTuple, seq);

  return make_tuple(make_kernel<KPol>(camp::get<Is>(lowBoundaryIterSpaceDecomp), camp::get<0>(knl.bodies))...);
}
  

// returns the kernels that execute the lower boundary iterations for a loopchain
// Within the loop chain, there are iteration shared by all the loops in the chain.
// This function generates the kernels that execute the iterations that will be 
// executed before the shared iterations
template <typename...KernelTypes, camp::idx_t...Is>
auto low_boundary_knls(camp::tuple<KernelTypes...> knlTuple, camp::idx_seq<Is...> seq) {
  
  auto segmentTuples = make_tuple(camp::get<Is>(knlTuple).segments...);
  auto sharedIterSpaceTuple = intersect_segment_tuples(segmentTuples);

  return tuple_cat(low_bondary_knls_for_knl(camp::get<Is>(knlTuple), 
                                            sharedIterSpaceTuple, 
                                            idx_seq_for(sharedIterSpaceTuple))...);
}










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
  return fuse(knlTuple, idx_seq_for(knlTuple);
}

template <typename...KernelTypes>
auto fuse(KernelTypes... knls) {
  return fuse(make_tuple(knls));
}

} //namespace RAJA

#endif
