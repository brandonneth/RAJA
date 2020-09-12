// Contains common code for kernel transformations

#ifndef RAJA_LCCOMMON_HPP
#define RAJA_LCCOMMON_HPP

#include "RAJA/config.hpp"


namespace RAJA
{

// Function Declarations

template <typename...SegmentTupleTypes>
auto intersect_segment_tuples(camp::tuple<SegmentTupleTypes...> segmentTuples);

template <typename...SegmentTupleTypes>
auto intersect_segment_tuples(SegmentTupleTypes... segmentTuples);



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
  auto groupedByDim =  tuple_zip(segmentTuples);
  return make_tuple(intersect_segments(camp::get<NumDims>(groupedByDim))...);
}



template <typename...SegmentTupleTypes, camp::idx_t...Is>
auto intersect_segment_tuples(camp::tuple<SegmentTupleTypes...> segmentTuples, camp::idx_seq<Is...> seq){
  return intersect_segment_tuples(segmentTuples, seq, idx_seq_for(camp::get<0>(segmentTuples)));
}
// For an arbitrary number of segment tuples of the same dimensionality,
// returns a segment tuple for the intersection of the segment tuples.
// This is the shared iteration space of the kernels with the provided segment tuple
template <typename...SegmentTupleTypes>
auto intersect_segment_tuples(camp::tuple<SegmentTupleTypes...> segmentTuples) {
  return intersect_segment_tuples(segmentTuples, idx_seq_for(segmentTuples));
}


template <typename...SegmentTupleTypes>
auto intersect_segment_tuples(SegmentTupleTypes... segmentTuples) {
  return intersect_segment_tuples(make_tuple(segmentTuples...));
}

// Common ISL Code Declarations


template <typename KernelType>
isl_union_set * knl_iterspace(isl_ctx * ctx, KernelType knl, camp::idx_t id);

template <typename KernelType>
isl_union_map * knl_read_relation(isl_ctx * ctx, KernelType knl, camp::idx_t id);

template <typename KernelType>
isl_union_map * knl_write_relation(isl_ctx * cxt, KernelType knl, camp::idx_t id);

template <typename KernelType1, typename KernelType2>
isl_union_map * flow_dep_relation(isl_ctx * ctx,
                                  KernelType1 knl1, KernelType2 knl2,
                                  camp::idx_t id1, camp::idx_t id2);

template <typename KernelType1, typename KernelType2>
isl_union_map * anti_dep_relation(isl_ctx * ctx,
                                  KernelType1 knl1, KernelType2 knl2,
                                  camp::idx_t id1, camp::idx_t id2);

template <typename KernelType1, typename KernelType2>
isl_union_map * output_dep_relation(isl_ctx * ctx,
                                  KernelType1 knl1, KernelType2 knl2,
                                  camp::idx_t id1, camp::idx_t id2);

template <typename KernelType1, typename KernelType2>
isl_union_map * data_dep_relation(isl_ctx * ctx,
                                  KernelType1 knl1, KernelType2 knl2,
                                  camp::idx_t id1, camp::idx_t id2);


// Common ISL Code Definitions


template <typename SegmentType>
std::string constraints_from_segment(SegmentType segment, camp::idx_t idx) {
  auto low = *segment.begin();
  auto high = *segment.end();

  std::stringstream s;
  s << "i" << idx << " >= " << low << " and " << "i" << idx << " < " << high;
  std::string str = s.str();

  return str;
}


template <typename...SegmentTypes, camp::idx_t I, camp::idx_t...Is>
auto segment_tuple_to_constraints(camp::tuple<SegmentTypes...> segmentTuple, camp::idx_seq<I, Is...> seq) {
   std::string rest;
  if constexpr (sizeof...(Is) > 0) {
    rest = segment_tuple_to_constraints<Is...>(segmentTuple, camp::idx_seq<Is...>{});
    return constraints_from_segment(camp::get<I>(segmentTuple), I) + " and " + rest;
  } else {
    return constraints_from_segment(camp::get<I>(segmentTuple));
  }
}



//returns the string for the codomain vector with numElements dimension. For example,
// 2 return "[i0,i1]" and 4 returns "[i0,i1,i2,i3]"
template <typename IdxType>
std::string codomain_vector(IdxType numElements, camp::idx_t loopNum) {

  std::stringstream vec;
  vec << "L" << loopNum << "[";

  for(int i = 0; i < numElements; i++) {
    vec << "i" << i;
    if(i < numElements - 1) {
      vec << ",";
    }
  }

  vec << "]";

  return vec.str();
}


template <typename KernelType>
isl_union_set * knl_iterspace(isl_ctx * ctx, KernelType knl, camp::idx_t id) {
  auto segments = knl.segments;

  auto bounds = segment_tuple_to_constraints(segments, idx_seq_for(segments));

  auto codomainString = codomain_vector(knl.numArgs, id);
 
  auto iterspaceString "{ " + codomainString + " : " + bounds + " }";

  isl_union_set * iterspace = isl_union_set_read_from_str(ctx, iterspaceString.c_str());

  return iterspace;
}


template <typename KernelType>
isl_union_map * read_relation(isl_ctx * ctx, KernelType knl, camp::idx_t id);

template <typename KernelType>
isl_union_map * write_relation(isl_ctx * cxt, KernelType knl, camp::idx_t id);

template <typename KernelType1, typename KernelType2>
isl_union_map * flow_dep_relation(isl_ctx * ctx,
                                  KernelType1 knl1, KernelType2 knl2,
                                  camp::idx_t id1, camp::idx_t id2);

template <typename KernelType1, typename KernelType2>
isl_union_map * anti_dep_relation(isl_ctx * ctx,
                                  KernelType1 knl1, KernelType2 knl2,
                                  camp::idx_t id1, camp::idx_t id2);

template <typename KernelType1, typename KernelType2>
isl_union_map * output_dep_relation(isl_ctx * ctx,
                                  KernelType1 knl1, KernelType2 knl2,
                                  camp::idx_t id1, camp::idx_t id2);

template <typename KernelType1, typename KernelType2>
isl_union_map * data_dep_relation(isl_ctx * ctx,
                                  KernelType1 knl1, KernelType2 knl2,
                                  camp::idx_t id1, camp::idx_t id2);



template <typename KernelType>

} //namespace RAJA
#endif

