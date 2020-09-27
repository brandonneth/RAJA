// Contains code for the overlapped tiling transformations

#ifndef RAJA_LCOVERLAPPEDTILE_HPP
#define RAJA_LCOVERLAPPEDTILE_HPP

#include "RAJA/config.hpp"
#include "RAJA/loopchain/KernelWrapper.hpp"
#include "RAJA/loopchain/Chain.hpp"
#include "RAJA/loopchain/transformations/Common.hpp"
namespace RAJA {


//creates the constraints for the overlap amounts between two kernels
//the overlap amount constraint is that for a particular dimension,
//the overlap amonut for knl1 has to be at least the amount for knl2
// plus the maximum length of the dependences between the kernels in
// that dimension. knl1_dim >= knl2_dim + maxDep(knl1, knl2, dim)
template <typename Knl1, typename Knl2>
std::string overlap_constraint(Knl1 knl1, Knl2 knl2, camp::idx_t id1, camp::idx_t id2) {
  assert(knl1.numArgs == knl2.numArgs);
  constexpr int numDims = knl1.numArgs;

  isl_ctx * ctx = isl_ctx_alloc();

  isl_union_map * dependences = data_dep_relation(ctx, knl1, knl2, id1, id2);
  
  isl_union_set * minInput = isl_union_set_lexmin(knl_iterspace(ctx, knl1, id1));

  isl_union_set * image = isl_union_set_appyl(isl_union_set_copy(minInput), dependences);

  if(isl_union_set_is_empty(image)) {
    return "true";
  }
  
  int * maxOut = new int[numDims + 1];
  int * inputValues = new int[numDims + 1];
  
  maxOut[0] = numDims;
  inputValues[0] = numDims;
  for(int i = 1; i < numDims; i++) {
    maxOut[i] = std::numeric_limits<int>::min();
    inputValues[i] = std::numeric_limits<int>::min();
  }

  // updates the maximum values of the dependence distances for each 
  // dimension of the dependences between the two kernels
  auto update_max_values = [] (isl_point * p, void * _currMax) {
    int * currMax = (int *) (currMax_);
    int numDims = currMax[0];
    for(int i = 0; i < numDims; i++) {
      isl_val * pointVal = isl_point_get_coordinate_val(p, isl_dim_set, i);
      int value = isl_val_get_num_si(pointVal);
      if(currMax[i+1] < value) {
        currMax[i+1] = value;
      }
    }
    return isl_stat_ok;
  };

  //use the func to fill out inputValues, which should only have one element
  isl_union_set_foreach_point(isl_union_set_copy(minInput), update_max_values, inputValues);
  // then use it to get the max points in the output set
  isl_union_set_foreach_point(isl_union_set_copy(image), update_max_values, maxOut);

  int * maxDependenceDistances = new int[numDims];
  for(int i = 0; i < numDims; i++) {
    maxDependenceDistance[i] = maxOut[i+1] - inputValues[i+1];
  }

  //last, build up the constraint string
  std::string constraints = "";
  for(int i = 0; i < numDims; i++) {
    std::string constraint = "";
    constraint += "O" + std::to_string(id1) + "_" + std::to_string(i);
    constraint += " >= ";
    constraint += "O" + std::to_string(id2) + "_" + std::to_string(i);
    constraint += "+";
    constraint += std::to_string(maxDependence[i]);
    constraints += constraint;
    if (i != numDims - 1) {constraints += " and ";}
  }
  
  return constraints;
} //overlap_constraint


template <camp::idx_t id1, camp::idx_t id2, camp::idx_t NumKnls, typename...KernelTypes>
auto overlap_constraints_helper(isl_ctx * ctx, camp::tuple<KernelTypes...> knlTuple) {
  if constexpr (KnlId1 > NumKnls) {
    std::string terminal = "END";
    return terminal;
  } else if constexpr (KnlId2 >= NumKnls) {
    return overlap_constraints_helper<id1+1, id1+2, NumKnls>(ctx, knlTuple);
  } else {
    auto constraintString = overlap_constraint(camp::get<id1>(knlTuple), camp::get<id2>(knlTuple), KnlId1, KnlId2);
    
    auto rest = overlap_constraints_helper<KnlId1, KnlId2+1, NumKnls>(ctx, knlTuple);

    if(rest == "END") {
      return constraintString;
    } else {
      return constraintString + " and " + rest;
    }
  }

}

template <typename...KernelTypes, camp::idx_t...Is>
auto overlap_constraints(isl_ctx * ctx, camp::tuple<KernelTypes...> knlTuple, camp::idx_seq<Is...> seq) {
  return overlap_constraints_helper<0,1,sizeof...(Is)>(ctx, knlTuple);
}

//performs the ISL calculation to determine the size of the overlaps for each tile.
//first, the strings for the space of all overlaps is created, then contrained by 
// a requirement of non-negativity and that the overlap amount between two loops 
// is at least as large as the biggest dependenec between them at that dimension
template <typename...KernelTypes, camp::idx_t...Is>
auto overlap_amount_tuples(camp::tuple<KernelTypes...> knlTuple, camp::idx_seq<Is...> seq) {
  
  constexpr auto numKernels = sizeof...(Is);
  assert(numKernels > 0);
  auto constexpr numDims = camp::get<0>(knlTuple).numArgs;

  isl_ctx * ctx = isl_ctx_alloc();

  std::string overlapSet = "[";
  std::string nonNegConstraint = "";
    for(int loopNum = 1; loopNum <= numKernels; loopNum++) {
    for(int dim = 0; dim < numDims; dim++) {
      std::string overlapString = "O" + std::to_string(loopNum) + "_" + std::to_string(dim);
      overlapSet += overlapString;

      nonNegConstraint += overlapString + " >= 0 ";

      if(dim != numDims - 1) {
       overlapSet += ",";
       nonNegConstraint += " and ";
      }
    }//dim
    if(loopNum != numKernels) {
      overlapSet += ",";
      nonNegConstraint += " and ";
    }
  }//loopNum

  overlapSet += "]";

  std::string legalOverlapsString = "{" + overlapSet + ": " + nonNegConstraint + " and " + overlap_constraints(ctx, knlTuple, seq) + "}";

  isl_union_set * legalOverlaps = isl_union_set_read_from_str(ctx, legalOverlapsString.c_str());

  isl_point * legalOverlap = isl_union_set_sample_point(isl_union_set_lexmin(legalOverlaps));

  auto point_to_vector = [](isl_point * p, int numDims) {
    std::vector<int> vals = {};

    for(int i = 0; i < numDims; i++) {
      isl_val * pointVal = isl_point_get_coordinate_val(p, isl_dim_set, i);
      int value = isl_val_get_num_si(pointVal);
      vals.push_back(value);
    }
    return vals;
  };

  auto overlapAmounts = point_to_vector(legalOverlap, numKernels * numDims);
  return shift_vector_to_shift_tuples<numKernels,numDims>(overlapAmounts);

}//overlap_amount_tuples



template <typename...KernelTypes, typename OverlapTupleType, camp::idx_t...Is>
auto overlapped_tile_no_fuse_executor(camp::tuple<KernelTypes...> knlTuple, 
                                      OverlapTupleType overlaps, 
                                      camp::idx_seq<Is...> seq) {


  auto entireRangeSegment = fused_segment_tuple(knlTuple, seq);

  auto f = [=](auto tileRangeSegments) {

     auto overlapTileSegmentTuples = make_tuple((overlap_segment_tuple(tileRangeSegments, entireRangeSegment, camp::get<Is>(overlaps)))...);
     auto newKnlTuple = make_tuple((change_segment_tuple(camp::get<Is>(knlTuple), camp::get<Is>(overlapTileSegmentTuples)))...);
     auto overlappedTileChain = chain(camp::get<Is>(newKnlTuple)...);
     overlappedTileChain();
  };

  return f;
}
auto overlapped_tile_no_fuse_executor(auto knlTuple, auto overlaps) {
  return overlapped_tile_no_fuse_executor(knlTuple, overlaps, idx_seq_for(knlTuple));
} //overlap_executor


//creates the kernel which tiles the overlapping region of the kernels in the tuple
template <camp::idx_t TileSize, camp::idx_t...Is>
auto overlap_tile_no_fuse_kernel(auto knlTuple, auto overlapAmountTuples, camp::idx_seq<Is...> knlSeq) {

  auto constexpr numDims = camp::get<0>(knlTuple).numArgs;
  using InnerPol = decltype(overlapped_tile_policy<TileSize,0,numDims>());
  using KPol = KernelPolicy<InnerPol>;

  auto segmentTuple = fused_segment_tuple(knlTuple, knlSeq);
  auto tiledExecutor = overlapped_tile_no_fuse_executor(knlTuple, overlapAmountTuples);

  return make_kernel<KPol>(segmentTuple, tiledExecutor);
}//overlap_kernel


//performs overlapped tiling without fusing the loops inside each tile
template <camp::idx_t TileSize, typename...KernelTypes, camp::idx_t...Is>
auto overlapped_tile_no_fuse(camp::tuple<KernelTypes...> knlTuple, camp::idx_seq<Is...> seq) {

  auto shiftAmountTuples = shift_amount_tuples(knlTuple, seq);

  auto shiftedKnlTuple = make_tuple(shift(camp::get<Is>(knlTuple), camp::get<Is>(shiftAmountTuples))...);

  auto preTilingKnls = low_boundary_knls(shiftedKnlTuple, seq);
  auto postTilingKnls = high_boundary_knls(shiftedKnlTuple, seq);


  auto overlapAmountTuples = overlap_amount_tuples(shiftedKnlTuple, seq);

  auto overlappedKernel = overlap_tile_no_fuse_kernel<TileSize, Is...>(shiftedKnlTuple, overlapAmountTuples, seq);

  auto allKernelTuple = tuple_cat(preTilingKnls, make_tuple(overlappedKernel), postTilingKnls);

  return grouped_kernels(allKernelTuple);
}//overlapped_tile_no_fuse



template <camp::idx_t TileSize, typename...KernelTypes>
auto overlapped_tile_no_fuse(KernelTypes... knls) {
  auto knlTuple = make_tuple(knls...);
  auto seq = idx_seq_for(knlTuple);

  return overlapped_tile_no_fuse<TileSize>(knlTuple, idx_seq_for(knlTuple));
}


template <typename...KernelTypes>
auto overlapped_tile_no_fuse(KernelTypes...knls) {
  constexpr auto TileSize = 16;
  auto knlTuple = make_tuple(knls...);
  return overlapped_tile_no_fuse<TileSize>(knlTuple, idx_seq_for(knlTuple));
} //overlapped_tile_no_fuse


} //namespace RAJA

#endif
