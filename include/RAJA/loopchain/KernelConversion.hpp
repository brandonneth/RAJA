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
std::string range_vector(int numElements) {

  std::stringstream vec;
  vec << "[";

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
template <typename KernelPol,typename SegmentTuple,typename...Bodies>
isl_union_set * iterspace_from_knl(const KernelWrapper<KernelPol,SegmentTuple,Bodies...> & knl) {

  isl_ctx* ctx;
  ctx = isl_ctx_alloc();
  isl_printer* p = isl_printer_to_file(ctx, stdout);

  auto segments = knl.segments;

  auto segmentSeq = idx_seq_for(segments);

  auto bounds = bounds_from_segment(segments, segmentSeq);

  auto rangeVector = range_vector(knl.numArgs);  

  auto iterspaceString = "{ " + rangeVector + " : " + bounds + " }";
 
  isl_union_set * iterspace = isl_union_set_read_from_str(ctx, iterspaceString.c_str());

  return iterspace; 

} // iterspace_from_knl


//Returns the data space accessed by a knl 
//











template <typename KernelPol,typename SegmentTuple,typename...Bodies>
isl_union_set * dataspace_from_knl(KernelWrapper<KernelPol,SegmentTuple,Bodies...> knl) {

  isl_ctx* ctx;
  ctx = isl_ctx_alloc();
  isl_printer* p = isl_printer_to_file(ctx, stdout);

  auto accesses = knl.execute_symbolically();

  for(auto access : accesses) {
    
  }
  read_relation_from_knl(knl,ctx);
  return NULL; 

} // dataspace_from_knl

template <typename KernelPol,typename SegmentTuple,typename...Bodies>
isl_union_map * read_relation_from_knl(KernelWrapper<KernelPol,SegmentTuple,Bodies...> knl, isl_ctx * ctx) {

  isl_printer* p = isl_printer_to_file(ctx, stdout);

  auto accesses = knl.execute_symbolically();
 
  std::stringstream readString;

  readString << "{";

  for(auto access : accesses) {
    if(access.isRead) {

      auto name = get_array_name(access);
      auto accessString = access.access_string();
      readString << "\t" << range_vector(knl.numArgs) << " -> " << name << "[" << accessString << "];";
    }
  }

  readString << "}";

  std::cout << "Read string for knl: " << readString.str() << "\n";
  
  isl_union_map * readRelation = isl_union_map_read_from_str(ctx, readString.str().c_str());

  p = isl_printer_print_union_map(p, readRelation);

  return NULL;

  
} // read_relation_from_knl






} //namespace RAJA
#endif
