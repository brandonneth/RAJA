// Contains functions that extract isl polyhedral sets from kernel objects.

#ifndef RAJA_LC_knlconv_HPP
#define RAJA_LC_knlconv_HPP

#include "RAJA/util/all-isl.h"
#include <iostream>
#include <sstream>

#include "RAJA/loopchain/KernelWrapper.hpp"
#include "RAJA/loopchain/Utils.hpp"
namespace RAJA {

std::string get_array_name(SymAccess a);







template <camp::idx_t Idx>
std::string bound_from_segment(auto segment) {
  
  auto low = *segment.begin();
  auto high = *segment.end();

  std::stringstream s;
  s << "i" << Idx << " >= " << low << " and " << "i" << Idx << " < " << high;
  std::string str = s.str();

  return str;
}

template <camp::idx_t I>
std::string bounds_from_segment(auto segmentTuple, camp::idx_seq<I>) {

  return bound_from_segment<I>(camp::get<I>(segmentTuple));
  
}

template <camp::idx_t I, camp::idx_t...Is>
std::string bounds_from_segment(auto segmentTuple, camp::idx_seq<I, Is...>) {
  std::string rest;
  if constexpr (sizeof...(Is) > 0) {
    rest = bounds_from_segment<Is...>(segmentTuple, camp::idx_seq<Is...>{});
    return bound_from_segment<I>(camp::get<I>(segmentTuple)) + " and " + rest;
  } else {
    return bound_from_segment<I>(camp::get<I>(segmentTuple));
  }
}


//returns the string for the range vector with numElements dimension. For example,
// 2 return "[i0,i1]" and 4 returns "[i0,i1,i2,i3]"
std::string range_vector(int numElements, int loopNum) {

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
} //range_vector

//Returns the iteration space polyhedron of knl. 
//Does so without the loop number in the chain, 
//so that must be prepended if necessary.
template <camp::idx_t LoopNum>
isl_union_set * iterspace_from_knl(isl_ctx * ctx, auto knl) {


  auto segments = knl.segments;

  auto segmentSeq = idx_seq_for(segments);

  auto bounds = bounds_from_segment(segments, segmentSeq);

  auto rangeVector = range_vector(knl.numArgs, LoopNum);  

  auto iterspaceString = "{ " + rangeVector + " : " + bounds + " }";
 
  isl_union_set * iterspace = isl_union_set_read_from_str(ctx, iterspaceString.c_str());

  return iterspace; 

} // iterspace_from_knl

template <camp::idx_t LoopNum>
isl_union_map * read_relation_from_knl(isl_ctx * ctx, auto knl) {

  isl_printer* p = isl_printer_to_file(ctx, stdout);

  auto accesses = knl.execute_symbolically();
 
  std::stringstream readString;

  readString << "{";

  for(auto access : accesses) {
    if(access.isRead) {

      auto name = get_array_name(access);
      auto accessString = access.access_string();
      readString << "\t" << range_vector(knl.numArgs, LoopNum) << " -> " << name << "[" << accessString << "];";
    }
  }

  readString << "}";

  
  isl_union_map * readRelation = isl_union_map_read_from_str(ctx, readString.str().c_str());

  isl_union_map * readRelationForLoop = isl_union_map_intersect_domain(readRelation, iterspace_from_knl<LoopNum>(ctx, knl));

  return readRelationForLoop;

  
} // read_relation_from_knl

template <camp::idx_t LoopNum>
isl_union_map * write_relation_from_knl(isl_ctx * ctx, auto knl) {

  isl_printer* p = isl_printer_to_file(ctx, stdout);

  auto accesses = knl.execute_symbolically();
 
  std::stringstream writeString;

  writeString << "{";

  for(auto access : accesses) {
    if(access.isWrite) {

      auto name = get_array_name(access);
      auto accessString = access.access_string();
      writeString << "\t" << range_vector(knl.numArgs, LoopNum) << " -> " << name << "[" << accessString << "];";
    }
  }

  writeString << "}";

  
  isl_union_map * writeRelation = isl_union_map_read_from_str(ctx, writeString.str().c_str());

  isl_union_map * writeRelationForLoop = isl_union_map_intersect_domain(writeRelation, iterspace_from_knl<LoopNum>(ctx, knl));

  return writeRelationForLoop;

  
} // read_relation_from_knl


//returns a relation between the iteration spaces of the two knls mapping instances of knl1 to 
//instances of knl2 where knl1 instance writes and knl2 instance reads the same location
template <camp::idx_t LoopNum1, camp::idx_t LoopNum2>
isl_union_map * flow_dep_relation_from_knls(isl_ctx * ctx,
                                            auto knl1, 
                                            auto knl2) {

  isl_union_map * knl1Writes = write_relation_from_knl<LoopNum1>(ctx, knl1);
  isl_union_map * knl2Reads = read_relation_from_knl<LoopNum2>(ctx, knl2);

  //We want a map from I1 to I2 where i1,i2 iff knl1WriteRelation(i1) and knl2ReadRelation(i2) have shared elements
 
  isl_union_map * readsInverse = isl_union_map_reverse(knl2Reads);

  //We apply the range of knl1Writes to the inverse of knl2Reads. 
  //The range of knl1Writes is the data accessed, which is the domain of the knl2Reads inverse
  //This basically feeds the data locations into the inversed knl2Reads, 
  //which gives out the iterations of knl2 that made the accesses
  isl_union_map * readsAfterWrites = isl_union_map_apply_range(knl1Writes, readsInverse);
  return readsAfterWrites;
} // flow_dep_relation_from_knls


//returns a relation between the iteration spaces of the two knls mapping instances of knl1 to 
//instances of knl2 where knl1 instance reads and knl2 instance writes the same location
template <camp::idx_t LoopNum1, camp::idx_t LoopNum2>
isl_union_map * anti_dep_relation_from_knls(isl_ctx * ctx,
                                            auto knl1, 
                                            auto knl2) {

  isl_union_map * knl1Reads = read_relation_from_knl<LoopNum1>(ctx, knl1);
  isl_union_map * knl2Writes = write_relation_from_knl<LoopNum2>(ctx, knl2);

 
  isl_union_map * writesInverse = isl_union_map_reverse(knl2Writes);
  isl_union_map * writesAfterReads = isl_union_map_apply_range(knl1Reads, writesInverse);
  return writesAfterReads;
} // anti_dep_relation_from_knls

//returns a relation between the iteration spaces of the two knls mapping instances of knl1 to 
//instances of knl2 where knl1 instance writes and knl2 instance writes the same location
template <camp::idx_t LoopNum1, camp::idx_t LoopNum2>
isl_union_map * output_dep_relation_from_knls(isl_ctx * ctx,
                                            auto knl1, 
                                            auto knl2) {

  isl_union_map * knl1Writes = write_relation_from_knl<LoopNum1>(ctx, knl1);
  isl_union_map * knl2Writes = write_relation_from_knl<LoopNum2>(ctx, knl2);

 
  isl_union_map * writesInverse = isl_union_map_reverse(knl2Writes);
  isl_union_map * writesAfterWrites = isl_union_map_apply_range(knl1Writes, writesInverse);
  return writesAfterWrites;
} // output_dep_relation_from_knls

//returns a relation between the iteration spaces of the two knls mapping instances of knl1 to instances of knl2 where there is a data dependence of some sort between the two instances.
template <camp::idx_t LoopNum1, camp::idx_t LoopNum2>
isl_union_map * data_dep_relation_from_knls(isl_ctx * ctx,
                                            auto knl1, 
                                            auto knl2) {

  isl_union_map * raw = flow_dep_relation_from_knls<LoopNum1,LoopNum2>(ctx, knl1, knl2);
  isl_union_map * war = anti_dep_relation_from_knls<LoopNum1,LoopNum2>(ctx, knl1, knl2);
  isl_union_map * waw = output_dep_relation_from_knls<LoopNum1,LoopNum2>(ctx, knl1, knl2);

  isl_union_map * partialUnion = isl_union_map_union(war, waw); 
  isl_union_map * depUnion = isl_union_map_union(raw, partialUnion);

  return depUnion;

} // data_dep_relation_from_knls


//returns the range string for normal execution of a loop with numdims that is loop loopNum in teh chain
std::string original_schedule_range(int numDims, int loopNum) {

  std::stringstream s;

  s << "[";
  s << loopNum << ",";

  for(int i = 0; i < numDims; i++) {
    s << "i" << i << ",0";
    if( i != numDims - 1) {
      s << ",";
    }
  }
  
  s << "]";

  return s.str();
} // original_schedule_domain

// Returns the schedule for the loop executed sequentially. For example,
// a double nested loop will give [i,j] -> [loopNum, i, 0, j, 0]
template <camp::idx_t LoopNum>
isl_union_map * original_schedule(isl_ctx * ctx, auto knl) {

  isl_union_set * iterspace = iterspace_from_knl<LoopNum>(ctx, knl);

  std::string scheduleRange = original_schedule_range(knl.numArgs, LoopNum);

  
  
  std::string scheduleRelation = "{" + range_vector(knl.numArgs, LoopNum) + " -> " + scheduleRange + "}";

  isl_union_map * scheduleNoBounds = isl_union_map_read_from_str(ctx, scheduleRelation.c_str());

  isl_printer * p = isl_printer_to_file(ctx, stdout);

  std::cout << "\nSchedule relation no bounds:";
  p = isl_printer_print_union_map(p, scheduleNoBounds);

  isl_union_map * scheduleWithBounds = isl_union_map_intersect_domain(scheduleNoBounds, iterspace);

  std::cout << "\nSchedule with bounds:";
  p = isl_printer_print_union_map(p, scheduleWithBounds);
   
  
  return scheduleWithBounds;
} //original_schedule



} //namespace RAJA
#endif
