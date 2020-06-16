// Contains the code for the transformation directives

#ifndef RAJA_LCTRANSFORMATIONS_HPP 
#define RAJA_LCTRANSFORMATIONS_HPP
#include "RAJA/config.hpp"


#include "RAJA/loopchain/Utils.hpp"
namespace RAJA
{



template <camp::idx_t LoopId, typename...ShiftAmounts>
struct Shift {
  
  camp::tuple<ShiftAmounts...> amounts;
  
  Shift(ShiftAmounts... _amounts) : amounts(_amounts...) {}
  
}; //Shift 

template <camp::idx_t LoopId, typename...ShiftAmounts>
auto shift(ShiftAmounts&&... amounts) {
  return Shift<LoopId, ShiftAmounts...>(std::forward<ShiftAmounts>(amounts)...);
}


template <camp::idx_t I >
auto shift_dimension(const auto & segmentTuple, auto shiftTuple) {
  
  auto dimension = camp::get<I>(segmentTuple);
  auto shiftAmount = camp::get<I>(shiftTuple);
  
  using RangeType = decltype(dimension);

  auto low = *dimension.begin();
  auto high = *dimension.end();
 
  return RangeType(low+shiftAmount, high+shiftAmount);
  
}

template <camp::idx_t... Is>
auto shift_dimensions(const auto & segmentTuple, auto shiftTuple, camp::idx_seq<Is...>) {
  
  auto shiftedIterationSpace = camp::make_tuple(shift_dimension<Is>(segmentTuple, shiftTuple)...);
  return shiftedIterationSpace; 
}

template <camp::idx_t... Is>
auto shift_bodies(auto bodies, auto shiftTuple, camp::idx_seq<Is...>) {
  return 0;
}

//hardcoding for 1 2 and 3 dimensional loops
template <camp::idx_t I1>
auto shift_body(const auto & body, auto shiftTuple, camp::idx_seq<I1>) {

  auto newBody = [=](auto newI) {
    body(newI-camp::get<0>(shiftTuple));
  };

  return newBody;
  
}

template<typename KernelPol, typename Segment, typename...Bodies>
auto shift_kernel(const KernelWrapper<KernelPol,Segment,Bodies...> & knl, auto shift) {
  
  auto iterationSpaceTuple = knl.segments;  
  
  auto predim1 = camp::get<0>(iterationSpaceTuple);
  auto shiftAmountTuple = shift.amounts;

  auto iterationSpaceSeq = camp::make_idx_seq_t<camp::tuple_size<decltype(iterationSpaceTuple)>::value>{};

  auto newIterationSpaceTuple = shift_dimensions(iterationSpaceTuple, shiftAmountTuple, iterationSpaceSeq);

  auto dim1 = camp::get<0>(newIterationSpaceTuple);


  auto shiftedBody = shift_body(camp::get<0>(knl.bodies), shiftAmountTuple, iterationSpaceSeq);
  return make_kernel<KernelPol>(newIterationSpaceTuple, shiftedBody);
   
}



template <camp::idx_t... LoopIds>
struct Fusion {

}; //Fusion

template<camp::idx_t... LoopIds>
auto fuse() {
  return Fusion<LoopIds...>();
}


template <camp::idx_t I>
auto calculate_overlap_segments_helper(auto seg1, auto seg2) {
  
  auto dim1 = camp::get<I>(seg1);
  auto dim2 = camp::get<I>(seg2);

  using RangeType = decltype(dim1);

  auto low = std::max(*dim1.begin(), *dim2.begin());

  auto high = std::min(*dim1.end(), *dim2.end());

 
  return RangeType(low,high); 
}

//takes two kernels and returns the segment tuple for their overlap
template <camp::idx_t... Is>
auto calculate_overlap_segments(auto seg1, auto seg2, camp::idx_seq<Is...>) {

  auto tuple = make_tuple((calculate_overlap_segments_helper<Is>(seg1, seg2))...);  

  return tuple;
}


template <camp::idx_t I>
auto pre_fuse_first_ith_dims_helper(auto originalSegments, auto overlapSegments) {

  auto fuseSegment = camp::get<I>(overlapSegments);
  auto originalSegment = camp::get<I>(originalSegments);

  using RangeType = decltype(fuseSegment);

  return RangeType(*fuseSegment.begin(), *originalSegment.end());

}

template <camp::idx_t... Is>
auto pre_fuse_first_ith_dims(auto originalSegments, auto overlapSegments, camp::idx_seq<Is...>) {
  
  return make_tuple((pre_fuse_first_ith_dims_helper<Is>(originalSegments, overlapSegments))...);  
}

template <camp::idx_t... Is>
auto pre_fuse_last_dims(auto originalSegments, camp::idx_seq<Is...>) {

  return make_tuple(camp::get<Is>(originalSegments)...);
}
//returns the segment tuple for the I-th pre-fuse kernel
template <camp::idx_t I, camp::idx_t... Is>
auto pre_fuse_kernels_segments(auto originalSegments, auto overlapSegments, camp::idx_seq<Is...>) {
  
  //up to the I-th dimension, it is fuse.start to original.end

  auto firstIthSeq = idx_seq_from_to<0,I>();
  auto firstIthDims = pre_fuse_first_ith_dims(originalSegments, overlapSegments, firstIthSeq);
  
  // for the I-th, its orig.start to fuse.start

  auto ithOriginal = camp::get<I>(originalSegments);
  auto ithFuse = camp::get<I>(overlapSegments);

  using RangeType = decltype(ithOriginal);

  auto ithDim = make_tuple(RangeType(*ithOriginal.begin(), *ithFuse.begin()));

  // for ith to end, its original.start to original.end

  auto endSeq = idx_seq_from_to<I,sizeof...(Is)-1>();
  auto endDims = pre_fuse_last_dims(originalSegments, endSeq);

  return camp::tuple_cat_pair(firstIthDims, camp::tuple_cat_pair(ithDim, endDims));
  
  
  
}

template <typename KernelPol, camp::idx_t ...Is>
auto kernel_helper(auto segmentTuple, auto bodiesTuple, camp::idx_seq<Is...>) {
  return make_kernel<KernelPol>(segmentTuple, camp::get<Is>(bodiesTuple)...);
}

template <typename KernelPol, typename Segments, typename...Bodies, camp::idx_t... Is>
auto pre_fuse_kernels(KernelWrapper<KernelPol,Segments,Bodies...> knl, auto overlapSegments, camp::idx_seq<Is...> seq) {

  auto bodies = knl.bodies;

  auto seg = knl.segments;

  auto segmentTuples = make_tuple(pre_fuse_kernels_segments<Is>(seg, overlapSegments, seq)...);
  
  return make_tuple((make_kernel<KernelPol>(camp::get<Is>(segmentTuples), camp::get<0>(bodies)))...);
   
} //pre_fuse_kernels

template<camp::idx_t I>
auto post_fuse_last_dims_helper(auto originalSegments, auto overlapSegments) {
  
  auto fuseSegment = camp::get<I>(overlapSegments);
  auto originalSegment = camp::get<I>(originalSegments);

  using RangeType = decltype(fuseSegment);

  return RangeType(*fuseSegment.begin(), *originalSegment.end());


}

template <camp::idx_t...Is>
auto post_fuse_last_dims(auto originalSegments, auto overlapSegments, camp::idx_seq<Is...>) {

  
  return make_tuple((post_fuse_last_dims_helper<Is>(originalSegments, overlapSegments))...);  
}

//returns the segment tuple for the I-th post-fuse kernel
template <camp::idx_t I, camp::idx_t... Is>
auto post_fuse_kernels_segments(auto originalSegments, auto overlapSegments, camp::idx_seq<Is...>) {
  
  //up to the I-th dimension, it is fuse.start to fuse.end

  auto firstIthDims = slice_tuple<0,I>(overlapSegments);


  // for ith, it is fuse.end to original.end
  auto ithOriginal = camp::get<I>(originalSegments);
  auto ithFuse = camp::get<I>(overlapSegments);

  using RangeType = decltype(ithOriginal);

  auto ithDim = make_tuple(RangeType(*ithFuse.end(), *ithOriginal.end()));

  // for ith to end, its fuse.start to original.end

  auto endSeq = idx_seq_from_to<I,sizeof...(Is)-1>();
  auto endDims = post_fuse_last_dims(originalSegments,overlapSegments, endSeq);

  
  auto postFuseSegments = camp::tuple_cat_pair(firstIthDims, camp::tuple_cat_pair(ithDim, endDims));
  
  auto size = camp::tuple_size<decltype(postFuseSegments)>().value;

   return postFuseSegments; 
}


template <typename KernelPol, typename Segments, typename...Bodies, camp::idx_t... Is>
auto post_fuse_kernels(KernelWrapper<KernelPol,Segments,Bodies...> knl, auto overlapSegments, camp::idx_seq<Is...> seq) {
  auto bodies = knl.bodies;

  auto seg = knl.segments;

  auto segmentTuples = make_tuple(post_fuse_kernels_segments<Is>(seg, overlapSegments, seq)...);

  return make_tuple(make_kernel<KernelPol>(camp::get<Is>(segmentTuples), camp::get<0>(bodies))...);


}


template <typename KernelPol, camp::idx_t I1>
auto fused_kernel(auto body1, auto body2, auto overlapSegment, camp::idx_seq<I1>) {

  auto lambda = [=](auto i) {
    body1(i);
    body2(i);
  };
  
  
  return make_kernel<KernelPol>(overlapSegment,lambda);
}

template <typename KernelPol, camp::idx_t I1, camp::idx_t I2>
auto fused_kernel(auto body1, auto body2, auto overlapSegment, camp::idx_seq<I1, I2>) {

  auto lambda = [=](auto i, auto j) {
    body1(i,j);
    body2(i,j);
  };
  
  return make_kernel<KernelPol>(overlapSegment,lambda);
}

template <typename KernelPol1, typename Segment1, typename...Bodies1, 
          typename KernelPol2, typename Segment2, typename...Bodies2>
auto fuse_kernels(const KernelWrapper<KernelPol1,Segment1,Bodies1...> & knl1, 
                  const KernelWrapper<KernelPol2,Segment2,Bodies2...> & knl2) {

  
  auto iterationSpaceTuple = knl1.segments;  
  auto iterationSpaceSeq = camp::make_idx_seq_t<camp::tuple_size<decltype(iterationSpaceTuple)>::value>{};
  auto overlapSegments = calculate_overlap_segments(knl1.segments, knl2.segments, iterationSpaceSeq);

  auto knl1PreKnls = pre_fuse_kernels(knl1, overlapSegments, iterationSpaceSeq);
  auto knl2PreKnls = pre_fuse_kernels(knl2, overlapSegments, iterationSpaceSeq);
  auto knl1PostKnls = post_fuse_kernels(knl1, overlapSegments, iterationSpaceSeq);
  auto knl2PostKnls = post_fuse_kernels(knl2, overlapSegments, iterationSpaceSeq);

  //TODO: Support multi-lambda kernels
  auto fusedKnl = fused_kernel<KernelPol1>(camp::get<0>(knl1.bodies), camp::get<0>(knl2.bodies), overlapSegments, iterationSpaceSeq);
  
  //return make_tuple(knl1,knl2);
  //return knl1PreKnls; 
  return tuple_cat(knl1PreKnls, knl2PreKnls, make_tuple(fusedKnl), knl1PostKnls, knl2PostKnls);
  
   

} //fuse_kernels



template <camp::idx_t LoopId, typename OverlapTupleType, typename TileSizeTupleType>
struct OverlappedTile {
  
  OverlapTupleType overlapAmounts;
  TileSizeTupleType tileSizes;

  OverlappedTile(auto _overlapAmounts, 
                 auto _tileSizes) : 
    overlapAmounts(_overlapAmounts), tileSizes(_tileSizes) {}

  //TODO: support different overlaps for different loops
  

}; //OverlappedTile

template < camp::idx_t LoopId, typename OverlapTupleType, typename TileSizeTupleType>
auto overlapped_tile(OverlapTupleType overlapAmounts, TileSizeTupleType tileSizes) {
  return OverlappedTile<LoopId, OverlapTupleType, TileSizeTupleType>(overlapAmounts, tileSizes);


} //overlapped_tile()

template <camp::idx_t LoopId, typename...TileSizeTypes>
auto tile(camp::tuple<TileSizeTypes...> tileSizes) {
  auto overlapAmounts = tuple_repeat<TileSizeTypes...>(0);

  return overlapped_tile<LoopId>(overlapAmounts, tileSizes);
}

template <camp::idx_t...Is>
std::vector<camp::idx_t> tuple_to_vector(auto tuple, camp::idx_seq<Is...>) {

  auto vector = std::vector<camp::idx_t>{camp::get<Is>(tuple)...};

  return vector;

}//tuple_to_vector

template <long int N, long int M, typename ExecPolicy>
auto extract_for_indices(RAJA::KernelPolicy<RAJA::statement::For<N,ExecPolicy,RAJA::statement::Lambda<M>>>) {
  return camp::idx_seq<N>{};
}
template <typename ExecPolicy, typename Statements, long int N>
auto extract_for_indices(RAJA::KernelPolicy<RAJA::statement::For<N,ExecPolicy,Statements>>) {

  using sub_policy_t = RAJA::KernelPolicy<Statements>;
  auto subPolicy = sub_policy_t();
  auto subIndices = extract_for_indices(subPolicy);

  auto thisIndex = camp::idx_seq<N>{};
  
  auto indices = idx_seq_cat(subIndices,thisIndex);

  return indices;
  
} //extract_for_indices;

template <typename Statements, camp::idx_t I>
auto apply_wrap(RAJA::KernelPolicy<Statements>, camp::idx_seq<I>) {
  using wrapped_policy_t = RAJA::KernelPolicy<RAJA::statement::OverlappedTile<I,0,RAJA::statement::tile_fixed<4>, RAJA::seq_exec, Statements>>;

  return wrapped_policy_t();
} //apply_wrap

template <typename Statements, camp::idx_t I, camp::idx_t...Is>
auto apply_wrap(RAJA::KernelPolicy<Statements>, camp::idx_seq<I,Is...>) {
  using wrapped_policy_t = RAJA::KernelPolicy<RAJA::statement::OverlappedTile<I,0,RAJA::statement::tile_fixed<4>, RAJA::seq_exec, Statements>>;

  return apply_wrap(wrapped_policy_t(), camp::idx_seq<Is...>{});
} //apply_wrap

template <typename ExecPolicy, typename Statements, long int N>
auto wrap_fors_with_overlapped_tile(RAJA::KernelPolicy<RAJA::statement::For<N, ExecPolicy, Statements>> knlPolicy) {
 
  auto forIndices = extract_for_indices(knlPolicy);
   
  auto wrapped = apply_wrap(knlPolicy, forIndices);

  return wrapped;
}//wrap_for_with_overlapped_tile


template <typename KernelPol, typename Segment, typename...Bodies>
auto overlapped_tile_kernel(KernelWrapper<KernelPol,Segment,Bodies...> knl, auto overlapAmounts, auto tileSizes) {

  auto overlappedTilePol = wrap_fors_with_overlapped_tile(KernelPol());

  using OverlappedTileKernelPolicy = decltype(overlappedTilePol);
  //TODO: support multi-lambda kernels
  auto newKnl = make_kernel<OverlappedTileKernelPolicy>(knl.segments, camp::get<0>(knl.bodies));

  newKnl.overlapAmounts = tuple_to_vector(overlapAmounts, idx_seq_for(overlapAmounts));
  newKnl.tileSizes = tuple_to_vector(tileSizes, idx_seq_for(tileSizes));
  return newKnl;
} // overlapped_tile_kernel


} //namespace RAJA

#endif
