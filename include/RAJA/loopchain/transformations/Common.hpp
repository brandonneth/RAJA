// Contains common code for kernel transformations

#ifndef RAJA_LCCOMMON_HPP
#define RAJA_LCCOMMON_HPP

#include "RAJA/config.hpp"


namespace RAJA
{

// Function Declarations


template <typename...SegmentTupleTypes>
auto intersect_segment_tuples(camp::tuple<SegmentTupleTypes...> segmentTuples);




// Function Definitions

// given a tuple of segments, returns one segment that contained by all the segments within the tuple
template <typename...SegmentTypes, camp::idx_t...Is>
auto intersect_segments(camp::tuple<SegmentTypes...> segments, camp::idx_seq<Is...> seq) {
  auto highestLow = max((*((camp::get<Is>(segments)).begin()))...);
  auto lowestHigh = min((*((camp::get<Is>(segments)).end()))...);

  return RangeSegment(highestLow, lowestHigh);
}

template <typename...SegmentTypes>
auto intersect_segments(camp::tuple<SegmentTypes...> segments) {
  return intersect_segments(segments, idx_seq_for(segments));
}

template <typename...SegmentTupleTypes, camp::idx_t...NumTuples, camp::idx_t...NumDims>
auto intersect_segment_tuples(camp::tuple<SegmentTupleTypes...> segmentTuples, 
                              camp::idx_seq<NumTuples...>,
                              camp::idx_seq<NumDims...>) {

  //for n in numDims
    //nth of each
  auto zeroethDims = make_tuple((camp::get<0>(camp::get<NumTuples>(segmentTuples))...));

  auto groupedByDim =  make_tuple((make_tuple((camp::get<NumDims>(camp::get<NumTuples>(segmentTuples))...)))...);

  return make_tuple(intersect_segments(camp::get<NumDims>(groupedByDim))...);
}



template <typename...SegmentTupleTypes, camp::idx_t...Is>
auto intersect_segment_tuples(camp::tuple<SegmentTupleTypes...> segmentTuples, camp::idx_seq<Is...> seq) {
  constexpr numDims = tuple_len(camp::get<0>(segmentTuples));

  return intersect_segment_tuples(segmentTuples, seq, idx_seq_for(camp::get<0>(segmentTuples)));
}
// For an arbitrary number of segment tuples of the same dimensionality,
// returns a segment tuple for the intersection of the segment tuples.
// This is the shared iteration space of the kernels with the provided segment tuple
template <typename...SegmentTupleTypes>
auto intersect_segment_tuples(camp::tuple<SegmentTupleTypes...> segmentTuples) {
  
  return intersect_segment_tuples(segmentTuples, idx_seq_for(segmentTuples));
}



} //namespace RAJA
#endif

