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

template <typename...SegmentTupleTypes, camp::idx_t...Is>
auto intersect_segment_tuples(camp::tuple<SegmentTupleTypes...> segmentTuples, camp::idx_seq<Is...> seq) {
  
  

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

