#ifndef RAJA_loopchain_kernel_HPP
#define RAJA_loopchain_kernel_HPP

#include "RAJA/config.hpp"

#include "RAJA/util/LoopChain.hpp"


namespace RAJA
{

/*
template <long int... Is, long int... Js>
auto combine_sequences(const std::integer_sequence<long int, Is...>&, const std::integer_sequence<long int, Js...>&) {

  return std::integer_sequence<long int, Is..., Js...>{};

}


template <long int N>
auto extract_for_indices(const RAJA::KernelPolicy<RAJA::statement::Lambda<N>>&) {
  return std::integer_sequence<long int>{};
} 

template <typename ExecPolicy, typename Statements, long int N >
auto extract_for_indices(RAJA::KernelPolicy<RAJA::statement::For<N,ExecPolicy,Statements>>) {

  using sub_policy_t = RAJA::KernelPolicy<Statements>;

  auto subPolicy = sub_policy_t();
  
  auto subIndices = extract_for_indices(subPolicy);

  auto thisIndex = std::integer_sequence<long int, N>{};

  auto combinedIndices = combine_sequences(subIndices, thisIndex);
  
  return combinedIndices;
}


template <long int... Is>
void print_sequence(std::integer_sequence<long int, Is...>) {

  ((std::cout << Is << " "),...);
 
  std::cout << "\n";

}

template <typename Statements>

auto apply_wrap(RAJA::KernelPolicy<Statements>, std::integer_sequence<long int>) {
  return RAJA::KernelPolicy<Statements>();
}

template <typename Statements, long int I, long int... Is>
auto apply_wrap(RAJA::KernelPolicy<Statements>, std::integer_sequence<long int, I, Is...>) {

  using wrapped_policy_t = RAJA::KernelPolicy<RAJA::statement::OverlappedTile<I, 0, RAJA::statement::tile_fixed<4>, RAJA::seq_exec, Statements>>;

  return apply_wrap(wrapped_policy_t(), std::integer_sequence<long int, Is...>{});

}

template <typename ExecPolicy, typename Statements, long int N>
auto wrap_for_with_overlapped_tile(RAJA::KernelPolicy<RAJA::statement::For<N, ExecPolicy, Statements>> p ) {

  using sub_policy_t = RAJA::KernelPolicy<Statements>;

  auto subPolicy = sub_policy_t();
 
  auto forIndices = extract_for_indices(RAJA::KernelPolicy<RAJA::statement::For<N,ExecPolicy,Statements>>());

  auto wrapped = apply_wrap(p, forIndices);
 
  print_sequence(forIndices);
  
  return wrapped;


  
}








 

template<std::size_t s>
SymIter makeSymIter() {

  std::string iterName = "i" + std::to_string(s);
  
  return SymIter(iterName);

}

template <typename PolicyType, typename SegmentTuple, typename... Bodies>
struct KernelW {
  std::vector<int> overlaps;
  const SegmentTuple &segments;
  std::tuple<Bodies...> bodies;
  static constexpr int numArgs = camp::tuple_size<SegmentTuple>::value;
  KernelW(SegmentTuple const & s, Bodies const &... b) : segments(s), bodies(b...) {
    overlaps = std::vector<int>();
  }

  KernelW(const KernelW & k) = default;
  KernelW(KernelW && k) = default;
  KernelW & operator=(const KernelW & k) = default;
  KernelW & operator=(KernelW && k) = default;

  template <std::size_t... Is>
  auto iterator_tuple(std::index_sequence<Is...>) {
    auto iters = std::tuple((makeSymIter<Is>())...);

    return iters;
  }

  std::vector<SymAccess> collect_accesses(SymIter iterator) {
    std::vector<SymAccess> accesses = std::vector<SymAccess>();
    for(long unsigned int i = 0; i < iterator.accesses->size(); i++) {
      accesses.push_back(iterator.accesses->at(i));
    }
    return accesses;
  }

  template <typename... Iterators>
  std::vector<SymAccess> collect_accesses(SymIter iterator, Iterators&&... iters) {

    std::vector<SymAccess> accesses = collect_accesses(std::forward<Iterators>(iters)...);

    std::vector<SymAccess> allAccesses = std::vector<SymAccess>();

    for(long unsigned int i = 0; i < accesses.size(); i++) {
      allAccesses.push_back(accesses.at(i));
    } 

    for(long unsigned int i = 0; i < iterator.accesses->size(); i++) {
      allAccesses.push_back(iterator.accesses->at(i));
    }

    return allAccesses;


  }
  
  template <std::size_t ...Is>
  std::vector<SymAccess> collect_accesses_from_tuple(auto iterators, std::index_sequence<Is...>) {
    return collect_accesses(std::get<Is>(iterators)...);

  }
  std::vector<SymAccess> execute_symbolically() {

    std::cout << "Symbolically executing kernel, there are " << numArgs << " arguments to the kernel\n";
    SymIter i = SymIter("i0");
    SymIter j = SymIter("i1");

    auto func = std::get<0>(bodies);
    
    auto iterators = iterator_tuple(std::make_index_sequence<numArgs>());
    
    std::apply(func, iterators);
    
    auto i0 = std::get<0>(iterators);

    auto accesses = collect_accesses_from_tuple(iterators, std::make_index_sequence<numArgs>()); 
    
  
    //for(long unsigned int i = 0; i < accesses1->size(); i++) {allAccesses.push_back(accesses1->at(i));}
   // for(long unsigned int i = 0; i < accesses2->size(); i++) {allAccesses.push_back(accesses2->at(i));}
    return accesses;
  }

  template <std::size_t... Is>
  void execute(std::index_sequence<Is...>) {


     util::PluginContext context{util::make_context<PolicyType>()};
    util::callPreLaunchPlugins(context);

    using segment_tuple_t = typename IterableWrapperTuple<camp::decay<SegmentTuple>>::type;

    auto params = RAJA::make_tuple();
    using param_tuple_t = camp::decay<decltype(params)>;

    if(overlaps.size() != 0){
      using loop_data_t = internal::LoopData<PolicyType, segment_tuple_t, param_tuple_t, camp::decay<Bodies>...>;
      loop_data_t loop_data(overlaps, make_wrapped_tuple(segments), params, std::get<Is>(bodies)...);
      
      RAJA_FORCEINLINE_RECURSIVE
      internal::execute_statement_list<PolicyType>(loop_data);

    } else {
  
      RAJA::kernel<PolicyType>(segments, std::get<Is>(bodies)...);
    }
  }
  
  void operator () () {
     auto seq = std::index_sequence_for<Bodies...>{};
    execute(seq);
  }

};


template <std::size_t... LoopNums > 
struct Fuse {

  Fuse()  {
      }


};



template <std::size_t... LoopNums> 
struct OverlappedTile {

  OverlappedTile(){
 }

};


template <typename PolicyType, typename SegmentTuple, typename... Bodies>
KernelW<PolicyType,SegmentTuple,Bodies...> makeKernel(SegmentTuple const & segment, Bodies const &... bodies) {

   return KernelW<PolicyType,SegmentTuple,Bodies...>((segment), (bodies)...);
}

template <typename PolicyType, typename SegmentTuple, typename... Bodies1, typename... Bodies2>
auto fuse_kernels(
   KernelW<PolicyType,SegmentTuple,Bodies1...> knl1,
   KernelW<PolicyType,SegmentTuple,Bodies2...> knl2)
{
  //std::cout << "fusing kernels\n";

  //static bool canFuse = 1;


  static auto newlambda = [=] (auto i, auto j) {
     std::get<0>(knl1.bodies)(i,j);
     std::get<0>(knl2.bodies)(i,j);
  };

  static auto newKernel = makeKernel<PolicyType>(knl1.segments,newlambda);

  //std::cout << "\nexecuting fused kernel from within fuseKernels\n";

  return newKernel;
}

template <typename ExecPol, typename Container, typename ...Bodies, typename...Kernels>
std::vector<std::vector<int>> amount_to_shift_kernels(KernelW<ExecPol,Container,Bodies...> knl1, Kernels&&... kernels) {

  std::cout << "Amount to shift kernels\n";

  auto accessLists = std::vector<std::vector<SymAccess>>();

  return amount_to_shift_kernels(knl1.numArgs, accessLists, knl1, std::forward<Kernels>(kernels)...);

}

template <typename ExecPol, typename Container, typename ...Bodies, typename...Kernels>
std::vector<std::vector<int>> amount_to_shift_kernels(int numArgs, const std::vector<std::vector<SymAccess>> &accessLists, KernelW<ExecPol,Container,Bodies...> knl,Kernels&&... kernels) {
   
  auto accesses = knl.execute_symbolically();
  
  auto newAccessLists = std::vector<std::vector<SymAccess>>();

  for(auto lst: accessLists) {
    newAccessLists.push_back(lst);
  }
 
  newAccessLists.push_back(accesses); 
  return amount_to_shift_kernels(numArgs, newAccessLists, std::forward<Kernels>(kernels)...);


}


isl_stat find_min_point_val(isl_point * pnt, void * user);
template <typename ...Kernels>
std::vector<std::vector<int>> amount_to_shift_kernels(int numArgs, const std::vector<std::vector<SymAccess>> &accessLists) {

  std::cerr << "Finished collecting access information, performing shift calculations\n";

  std::cerr << "Gathered access information for " << accessLists.size() << " kernels, each with " << numArgs << " dims\n";

  auto numKernels = accessLists.size();

  auto loopDims = numArgs;


  auto readString = read_string(accessLists, loopDims);
  auto writeString = write_string(accessLists, loopDims); 


  std::cerr << "read string: "<< readString << "\n";
  std::cerr << "write string: " << writeString << "\n";

  isl_ctx* ctx = isl_ctx_alloc();
  isl_printer * p = isl_printer_to_file(ctx, stdout);

  isl_union_map * reads = isl_union_map_read_from_str(ctx, readString.c_str());
  isl_union_map * writes = isl_union_map_read_from_str(ctx, writeString.c_str());

  isl_union_map * reads_inverse = isl_union_map_reverse(isl_union_map_copy(reads));
  isl_union_map * writes_inverse = isl_union_map_reverse(isl_union_map_copy(writes));

  isl_union_map * raw = isl_union_map_apply_range(isl_union_map_copy(writes), isl_union_map_copy(reads_inverse));
  isl_union_map * waw = isl_union_map_apply_range(isl_union_map_copy(writes), isl_union_map_copy(writes_inverse));
  isl_union_map * war = isl_union_map_apply_range(isl_union_map_copy(reads), isl_union_map_copy(writes_inverse));

  
  std::cerr << "\nwaw\n";
  p = isl_printer_print_union_map(p,waw);
  std::cerr << "\nraw\n";
  p = isl_printer_print_union_map(p,raw);
  std::cerr << "\nwar\n";
  p = isl_printer_print_union_map(p,war);



  isl_union_map * deps = isl_union_map_union(raw, waw);
  deps = isl_union_map_union(deps, war);
  
  std::cerr << "\nDependences\n";
  p = isl_printer_print_union_map(p,deps);

  std::string shiftSet = "[";
  std::string nonNegConstraint = "";
  for(long unsigned int i = 0; i < numKernels; i++) {
    for(int d = 0; d < loopDims; d++) {
      std::string shiftString = "S" + std::to_string(i) + "_" + std::to_string(d);
      
      shiftSet += shiftString;
      nonNegConstraint += shiftString + ">= 0 ";
      
      if(d != loopDims - 1) {
        shiftSet += ",";
      }
    }
    if(i != numKernels -1 ) {
      shiftSet += ",";
      nonNegConstraint += " and ";
    }
  }
  shiftSet += "]";
  std::cerr << "\nshift set" << shiftSet << "\n";
 
  std::string allPossibleShiftsString = "{" + shiftSet + ":" + nonNegConstraint + "}";
  std::cerr << "\nnonneg shift set string: " << allPossibleShiftsString << "\n";
  isl_union_set * allPossibleShifts = isl_union_set_read_from_str(ctx, allPossibleShiftsString.c_str());

  std::cerr << "\nallPossibleShifts\n";

  isl_set * constrainedShifts = isl_set_from_union_set(allPossibleShifts);

  
  p = isl_printer_print_set(p, constrainedShifts);
  
  for(long unsigned int i = 0; i < numKernels; i++) {
    for(long unsigned int j = i + 1; j < numKernels; j++) {
      std::string srcLoop = "{L" + std::to_string(i) + "[" ;
      std::string dstLoop = "{L" + std::to_string(j) + "[";
      for(int k = 0; k < loopDims; k++) {
        srcLoop += "0";
        dstLoop += "m" + std::to_string(k);
        if( k != loopDims -1) {
          srcLoop += ",";
          dstLoop += ",";
        }
      }
      srcLoop += "] }";
      dstLoop += "] }";
      std::cerr << "\nsource loop: " << srcLoop << "\n";
      std::cerr << "dest loop: " << dstLoop << "\n";
    
      isl_union_set * srcSet = isl_union_set_read_from_str(ctx, srcLoop.c_str());
      isl_union_set * dstSet = isl_union_set_read_from_str(ctx, dstLoop.c_str());
      std::cerr << "srcSet\n";
      p= isl_printer_print_union_set(p,srcSet);

      std::clog << "\ndstSet\n";
      p = isl_printer_print_union_set(p,dstSet);
      isl_union_set * depsFromSrc = isl_union_set_apply(srcSet, isl_union_map_copy(deps));
      
      std::clog << "\ndepsFromSrc\n";
      p = isl_printer_print_union_set(p,depsFromSrc);


      isl_union_set * src2dstUnionSet = isl_union_set_intersect(isl_union_set_copy(depsFromSrc), dstSet);
      std::clog << "\nsrc2dstUnionSet\n";
      p= isl_printer_print_union_set(p, src2dstUnionSet);
 
      int * minOffset = new int;
      *minOffset = 12345; 
      if(!isl_union_set_is_empty(src2dstUnionSet)) {
        isl_set * src2dst = isl_set_from_union_set(src2dstUnionSet);

        isl_set_foreach_point(src2dst, find_min_point_val, (void*) minOffset); 

      } 
      std::cout << "\nminoffset is: " << *minOffset << "\n";;

      if(*minOffset != 12345) {
        std::string shift1 = "S" + std::to_string(i) + "_" + "0";
        std::string shift2 = "S" + std::to_string(j) + "_" + "0";

        std::string constraint = shift2 + " - " + shift1 + ">=" + std::to_string(-1 * *minOffset);
        std::clog << "shift constraint: " << constraint << "\n";
        std::string constraintSetString = "{" + shiftSet + ": " + constraint + "}";
        isl_set * constraintSet = isl_set_read_from_str(ctx, constraintSetString.c_str());
        constrainedShifts = isl_set_intersect(isl_set_copy(constrainedShifts), constraintSet); 
      }

     
    }
  }


 isl_point * point = isl_set_sample_point(constrainedShifts);

  std::cout << "\nsample valid shift\n";
  p = isl_printer_print_point(p,point);

  std::vector<std::vector<int>> shiftVectors = std::vector<std::vector<int>>();

  for(long unsigned int i = 0; i < numKernels; i++) {
    std::vector<int> loopShifts = std::vector<int>();
    for(int d = 0; d < loopDims; d++) {
        
      int shiftAmountIndex = i * loopDims + d;
      isl_val * pointVal = isl_point_get_coordinate_val(point, isl_dim_set, shiftAmountIndex);
      
      int shiftAmount = isl_val_get_num_si(pointVal);
      loopShifts.push_back(shiftAmount);

     
    }
    shiftVectors.push_back(loopShifts);
  }

std::cout << "\n";
  return shiftVectors;
} //amount_to_shift_kernels
  





template <typename ExecPol, typename Container, typename...Bodies, typename... Kernels> 
auto overlapped_tile_kernels(KernelW<ExecPol,Container,Bodies...> knl, Kernels&&... knls) {

  static auto shiftAmounts = amount_to_shift_kernels(knl, std::forward<Kernels>(knls)...);

  for (auto shiftVector : shiftAmounts) {
    std::cout << "Shift amount: " ;
    for(auto shiftAmount : shiftVector) {
      std::cout << shiftAmount << ", ";
    }
    std::cout << "\n";
  }

  //auto shiftedKernels = apply_shift_kernels(shiftAmounts, knl, std::forward<Kernels>(knls)...);

  

}


*/
/*
template <typename...Kernels, std::size_t ...Is>
auto overlapped_tiling_2d_kernels_indexset(const std::tuple<Kernels...> & knls, const auto & overlapVectors, const auto tileDimensions, std::index_sequence<Is...>) {

  auto knl0 = std::get<0>(knls);
  auto iterSpaceTuple = knl0.segments;

  auto iterSpaceI = camp::get<0>(iterSpaceTuple);
  auto iterSpaceJ = camp::get<1>(iterSpaceTuple);
 
  auto iLength = (*iterSpaceI.end()) - (*iterSpaceI.begin());
  auto jLength = (*iterSpaceJ.end()) - (*iterSpaceJ.begin());

  auto iBegin = *iterSpaceI.begin();
  auto jBegin = *iterSpaceJ.begin();

  auto iTileCount = iLength / std::get<0>(tileDimensions);
  auto jTileCount = jLength / std::get<1>(tileDimensions);


  auto lambdas = std::tuple(std::get<0>(std::get<Is>(knls).bodies)...);
  
  for iTile in 0 to iTileCount
    iStart = iTileCount * iTileSize + iBegin
    iEnd = iStart + iTileSize
    for jTile in 0 to jTileCount
      jStart = jTileCount * jTileSize + jBegin
      jEnd = jStart + jTileSize
      for i in iStart to iEnd
         for j in jStart to jEnd
           blah
  
  
  RAJA::TypedIndexSet<RAJA::RangeSegment> iSet;
  RAJA::TypedIndexSet<RAJA::RangeSegment> jSet;


  for(int i = 0; i < iTileCount; i++) {
    int iStart = i * std::get<0>(tileDimensions) + iBegin;
    int iEnd = iStart + std::get<0>(tileDimensions);
    std::cout << "Adding i index set: " << iStart << ", " << iEnd << "\n";
    iSet.push_back(RAJA::RangeSegment(iStart,iEnd));
  }

  for(int j = 0; j < jTileCount; j++) {
    int jStart = j * std::get<1>(tileDimensions) + jBegin;
    int jEnd = jStart + std::get<1>(tileDimensions);
    std::cout << "Adding j index set: " << jStart << ", " << jEnd << "\n";
    jSet.push_back(RAJA::RangeSegment(jStart,jEnd));
  }

  using ISET_EXECPOL = RAJA::ExecPolicy<RAJA::seq_segit, RAJA::seq_exec>;
  using JSET_EXECPOL = RAJA::ExecPolicy<RAJA::seq_segit, RAJA::seq_exec>;

  using TILED_POL = RAJA::KernelPolicy<
    RAJA::statement::For<0, ISET_EXECPOL,
      RAJA::statement::For<1, JSET_EXECPOL,
        RAJA::statement::Lambda<0>
      >
    >
  >;

  //RAJA::forall<ISET_EXECPOL>(iSet, [=] (auto i) {std::cout << "forall with juts one index set, i = " << i << "\n";});
  std::cout << "Executing with index set\n";
  //RAJA::kernel<TILED_POL>(RAJA::make_tuple(iSet,jSet), std::get<0>(lambdas));
  
  return 0;


}


*/
/*
 * Creates the overlapped tiling, fused kernel for a tuple of 2d kernels. 
 * Does not yet handle the edges of the tile (not running overlapped when the tile is the edge tile)
 * 
 * Assumptions: 
 *  (1) The kernels have the same iteration spaces
 *  (2) The tile dimensions evenly split up the iteration space
 */
/*
template <typename... Kernels, std::size_t ...Is>
auto overlapped_tiling_2d_kernels(const std::tuple<Kernels...> & knls, const auto & overlapVectors, const auto tileDimensions, std::index_sequence<Is...>) {

  auto knl0 = std::get<0>(knls);
  auto iterSpaceTuple = knl0.segments;

  auto iterSpaceI = camp::get<0>(iterSpaceTuple);
  auto iterSpaceJ = camp::get<1>(iterSpaceTuple);
 
  auto iLength = (*iterSpaceI.end()) - (*iterSpaceI.begin());
  auto jLength = (*iterSpaceJ.end()) - (*iterSpaceJ.begin());

  



  std::cout << "Overlapped tiling for " << sizeof...(Kernels) << " 2 dimensional kernels\n";



  std::cout << "Tile size: \n";
  std::cout << "  i dimension: " << std::get<0>(tileDimensions) << "\n";
  std::cout << "  j dimension: " << std::get<1>(tileDimensions) << "\n";
 
  std::cout << "Iteration Space Dimensions:\n";
  std::cout << "  i dimension: " << iLength << "\n";
  std::cout << "  j dimension: " << jLength << "\n";


  std::cout << "Creating overlap for first knl\n";
  std::cout << "Amount of overlap for first knl:\n";
  
  auto overlap0 = std::get<0>(overlapVectors);
  auto overlap0i = std::get<0>(overlap0);
  auto overlap0j = std::get<1>(overlap0);

  std::cout << "i,j small square size: " << std::get<0>(overlap0) << "," << std::get<1>(overlap0) << "\n";
  std::cout << "i side rectangle: " << std::get<0>(tileDimensions)<< "," << overlap0j << "\n";
  std::cout << "j side rectangle: " << overlap0i << "," << std::get<1>(tileDimensions) << "\n";

  
  auto overlap0_rect1_low_i = 0 - std::get<0>(overlap0);  
  auto overlap0_rect1_low_j = 0 - std::get<1>(overlap0);

  auto overlap0_rect1_high_i = 0;
  auto overlap0_rect1_high_j = std::get<1>(tileDimensions);

  std::cout << "First Kernel's Overlap rectangle 1 bounds: " << overlap0_rect1_low_i << "," << overlap0_rect1_low_j << " to " << overlap0_rect1_high_i << "," << overlap0_rect1_high_j << "\n";
  
  auto overlap0_rect2_low_i = 0;
  auto overlap0_rect2_low_j = 0 - std::get<1>(overlap0);

  auto overlap0_rect2_high_i = std::get<0>(tileDimensions);
  auto overlap0_rect2_high_j = 0;


  std::cout << "First Kernel's Overlap rectangle 2 bounds: " << overlap0_rect2_low_i << "," << overlap0_rect2_low_j << " to " << overlap0_rect2_high_i << "," << overlap0_rect2_high_j << "\n";

  auto overlap0SegmentTuple0 = RAJA::make_tuple(
    RAJA::RangeSegment(overlap0_rect1_low_i, overlap0_rect1_high_i),
    RAJA::RangeSegment(overlap0_rect1_low_j, overlap0_rect1_high_j));

  auto overlap0SegmentTuple1 = RAJA::make_tuple(
    RAJA::RangeSegment(overlap0_rect2_low_i, overlap0_rect2_high_i),
    RAJA::RangeSegment(overlap0_rect2_low_j, overlap0_rect2_high_j));
    

  auto knl1 =  std::get<0>(knls);
  auto overlap_knl1 = [=](auto i, auto j, auto iTile, auto jTile) {
    std::get<0>(knl1.bodies)(iTile*std::get<0>(tileDimensions) + i,jTile * std::get<1>(tileDimensions) + j);
  };

  auto overlap_knls = std::tuple([=](auto i, auto j, auto iTile, auto jTile) {
    std::get<0>(std::get<Is>(knls).bodies) (iTile * std::get<0>(tileDimensions) + i, jTile * std::get<1>(tileDimensions) + j);
  }...);

  int numTiles_i = 5;
  int numTiles_j = 5;
 
  using OVERLAP1POLICY = RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::seq_exec,
      RAJA::statement::For<1, RAJA::seq_exec,
        RAJA::statement::For<2, RAJA::seq_exec,
          RAJA::statement::For<3, RAJA::seq_exec,
            RAJA::statement::Lambda<0, RAJA::statement::Segs<2,3,0,1>>
          >
        >,
        RAJA::statement::For<4, RAJA::seq_exec,
          RAJA::statement::For<5, RAJA::seq_exec,
            RAJA::statement::Lambda<0, RAJA::statement::Segs<4,5,0,1>>
          >
        >
      >
    >
  >;
*/
  /*
  std::cout << "Executing overlap for knl1\n";
  RAJA::kernel<OVERLAP1POLICY>(
    RAJA::make_tuple(
      RAJA::RangeSegment(0,numTiles_i),
      RAJA::RangeSegment(0,numTiles_j),
      RAJA::RangeSegment(overlap0_rect1_low_i, overlap0_rect1_high_i),
      RAJA::RangeSegment(overlap0_rect1_low_j, overlap0_rect1_high_j),
      RAJA::RangeSegment(overlap0_rect2_low_i, overlap0_rect2_high_i),
      RAJA::RangeSegment(overlap0_rect2_low_j, overlap0_rect2_high_j)),
    overlap_knl1);
  
*/

  /*



  auto overlaps_i = std::tuple(std::get<0>(std::get<Is>(overlapVectors))...);
  auto overlaps_j = std::tuple(std::get<1>(std::get<Is>(overlapVectors))...);
  
  auto overlaps_rect1_low_i = std::tuple((0 - std::get<Is>(overlaps_i))...);
  auto overlaps_rect1_low_j = std::tuple((0 - std::get<Is>(overlaps_j))...);

  auto overlaps_rect1_high_i = std::tuple((0 * std::get<Is>(overlaps_i))...); // need a tuple of zeros
  auto overlaps_rect1_high_j = std::tuple((std::get<1>(tileDimensions) + std::get<Is>(overlaps_j) - std::get<Is>(overlaps_j))...); // need a tuple of std::get<1>(tileDimensions)

  

  auto overlaps_rect2_low_i = std::tuple((0 * std::get<Is>(overlaps_i))...); // need a tuple of zeros
  auto overlaps_rect2_low_j = std::tuple((0 - std::get<Is>(overlaps_j))...);
  
  auto overlaps_rect2_high_i = std::tuple((std::get<0>(tileDimensions) + std::get<Is>(overlaps_i) - std::get<Is>(overlaps_i))...);

 
  auto overlaps_rect2_high_j = std::tuple((0 * std::get<Is>(overlaps_j))...); // need a tuple of zeros


  auto rect1_ibounds = std::tuple(RAJA::RangeSegment(std::get<Is>(overlaps_rect1_low_i), std::get<Is>(overlaps_rect1_high_i))...);
  auto rect1_jbounds = std::tuple(RAJA::RangeSegment(std::get<Is>(overlaps_rect1_low_j), std::get<Is>(overlaps_rect1_high_j))...);

  auto rect1_bounds = std::tuple(std::tuple(std::get<Is>(rect1_ibounds), std::get<Is>(rect1_jbounds))...);

  auto rect1_segments = RAJA::make_tuple(
    RAJA::RangeSegment(0,numTiles_i),
    RAJA::RangeSegment(0,numTiles_j),
    std::get<0>(std::get<Is>(rect1_bounds))...,
    std::get<1>(std::get<Is>(rect1_bounds))...);

  using RECT1POLICY = RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::seq_exec,
      RAJA::statement::For<1, RAJA::seq_exec,
        RAJA::statement::For<2+Is, RAJA::seq_exec,
          RAJA::statement::For<2+sizeof...(Kernels)+Is, RAJA::seq_exec,
            RAJA::statement::Lambda<Is, RAJA::statement::Segs<2+Is,2+sizeof...(Kernels)+Is,0,1>>
          >
        >...
      >
    >
  >;
       
  
  auto rect2_ibounds = std::tuple(RAJA::RangeSegment(std::get<Is>(overlaps_rect2_low_i), std::get<Is>(overlaps_rect2_high_i))...);

  auto rect2_jbounds = std::tuple(RAJA::RangeSegment(std::get<Is>(overlaps_rect2_low_j), std::get<Is>(overlaps_rect2_high_j))...);

  auto rect2_segments = RAJA::make_tuple(
    RAJA::RangeSegment(0,numTiles_i),
    RAJA::RangeSegment(0,numTiles_j),
    std::get<Is>(rect2_ibounds)...,
    std::get<Is>(rect2_jbounds)...);

  using RECT2POLICY = RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::seq_exec,
      RAJA::statement::For<1, RAJA::seq_exec,
        RAJA::statement::For<2+Is, RAJA::seq_exec,
          RAJA::statement::For<2+sizeof...(Kernels)+Is, RAJA::seq_exec,
            RAJA::statement::Lambda<Is, RAJA::statement::Segs<2+Is,2+sizeof...(Kernels)+Is,0,1>>
          >
        >...
      >
    >
  >;


  auto overlapping_tile_segments = RAJA::make_tuple(
    RAJA::RangeSegment(0,numTiles_i),
    RAJA::RangeSegment(0,numTiles_j),
    std::get<Is>(rect1_ibounds)...,
    std::get<Is>(rect1_jbounds)...,
    std::get<Is>(rect2_ibounds)...,
    std::get<Is>(rect2_jbounds)...,
    RAJA::RangeSegment(0,std::get<0>(tileDimensions)),
    RAJA::RangeSegment(0,std::get<1>(tileDimensions)));

  using OVERLAPPINGTILEPOLICY = RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::seq_exec,
      RAJA::statement::For<1, RAJA::seq_exec,
        RAJA::statement::For<2+Is, RAJA::seq_exec,
          RAJA::statement::For<2+sizeof...(Kernels)+Is, RAJA::seq_exec,
            RAJA::statement::Lambda<Is, RAJA::statement::Segs<2+Is,2+sizeof...(Kernels)+Is,0,1>>
          >
        >..., //overlap rectangle 1
        RAJA::statement::For<2+ 2 * sizeof...(Kernels)+Is, RAJA::seq_exec,
          RAJA::statement::For<2+3*sizeof...(Kernels)+Is, RAJA::seq_exec,
            RAJA::statement::Lambda<Is, RAJA::statement::Segs<2+2*sizeof...(Kernels)+Is,2+3*sizeof...(Kernels)+Is,0,1>>
          >
        >..., //overlap rectangle 2
        RAJA::statement::For<2+4*sizeof...(Kernels) + 0, RAJA::seq_exec,
          RAJA::statement::For<2+4*sizeof...(Kernels) + 1, RAJA::seq_exec,
            RAJA::statement::Lambda<Is, RAJA::statement::Segs<2+4*sizeof...(Kernels) + 0, 2+4*sizeof...(Kernels) + 1,0,1>>... // lambdas for tile
          >
        >
      > // jTile
    > // iTile
  >;

  auto overlappedKnl = RAJA::makeKernel<OVERLAPPINGTILEPOLICY>(overlapping_tile_segments, std::get<Is>(overlap_knls)...); 
  
  return overlappedKnl;
}//overlapped_tiling_2d_kernels


template <typename... Kernels, std::size_t ...Is, typename... ints> 
auto overlapped_tiling_1d_kernels(const std::tuple<Kernels...> & knlTuple, const std::tuple<ints...> overlapAmounts, const int tileSize, std::index_sequence<Is...>) {
  
  std::cout << "Creating overlapped tiling kernel with underlying tile size of " << tileSize << "\n";

  std::cout << "Extracting lambdas\n";

  auto lambdas = std::tuple((std::get<0>(std::get<Is>(knlTuple).bodies))...);
  

  std::cout << "Creating overlap lambdas\n";

  auto overlapLambdas = std::tuple(([=](auto i, auto tile) { (std::get<Is>(lambdas))(tile * tileSize + i);})...);

  std::cout << "Creating tile lambda\n";

  int tileLength = 8;
  int numTiles = tileLength / tileSize;
  auto tileSegment = RAJA::RangeSegment(0, numTiles);
 

  auto overlappedSegments = std::tuple(RAJA::RangeSegment(std::get<Is>(overlapAmounts) * -1, 0)...);
  
  auto tiledSegment = RAJA::RangeSegment(0,tileSize);

  using OVERLAPPEDPOL = RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::seq_exec, //tile loop
      RAJA::statement::For<Is+1, RAJA::seq_exec, //overlap1 loop
        RAJA::statement::Lambda<Is, RAJA::statement::Segs<Is+1,0>>>...,
      RAJA::statement::For<sizeof...(Kernels)+1, RAJA::seq_exec, // tiled loop part
        RAJA::statement::Lambda<Is+sizeof...(Kernels), RAJA::statement::Segs<sizeof...(Kernels)+1,0>>...
      >
    >
  >;

  auto totalTuple = RAJA::make_tuple(tileSegment,std::get<Is>(overlappedSegments)..., tiledSegment);

  auto tiledKnl = RAJA::makeKernel<OVERLAPPEDPOL>(totalTuple, std::get<Is>(overlapLambdas)...,std::get<Is>(overlapLambdas)...);

  std::cout << "Executing tiled kernel\n";
  tiledKnl();
  std::cout << "Returning tiled kernel\n"; 
  return tiledKnl; 
} //overlapped_tiling_1d_kernels



template <typename...KernelTemplate, std::size_t... indices>
auto shift_and_fuse_kernels(std::tuple<KernelTemplate...> knlTuple, std::index_sequence<indices...>) {

  auto amountToShift = amount_to_shift_kernels(std::get<indices>(knlTuple)...);

  
  
  return knlTuple;
}




template <typename...KernelTemplate, typename...Args> 
auto chain_kernels(KernelW<KernelTemplate...> knl, Args&&... args) {

  auto tuple = std::tuple(knl);

  return chain_kernels(tuple, std::forward<Args>(args)...);

}


template <typename... TupleTemplate, typename...KernelTemplate> 
auto chain_kernels(std::tuple<TupleTemplate...> tuple, KernelW<KernelTemplate...> knl) {

  std::cout << "chain_kernel case with no user provided directives.\n";

  std::cout << "TODO: Applying Default Trasnformation\n";

  static auto knlTuple = std::tuple(knl);
 
  static auto newTuple = std::tuple_cat(tuple, knlTuple);

  return newTuple; 

} 
template <typename... TupleTemplate, typename...KernelTemplate, typename...Args> 
auto chain_kernels(std::tuple<TupleTemplate...> tuple, KernelW<KernelTemplate...> knl, Args&&...args) {

  static auto knlTuple = std::tuple(knl);
 
  static auto newTuple = std::tuple_cat(tuple, knlTuple);

  return chain_kernels(newTuple, std::forward<Args>(args)...);

} 

template <std::size_t... LoopNums>
int is_in_sequence(std::size_t value) {

  std::vector<int> equalVector = std::vector<int>{(value == LoopNums)...};

  for(auto val : equalVector) {
    if(val){return 1;}
  }
  
  return 0;


}

//courtesy of https://stackoverflow.com/questions/53223910/how-to-access-n-th-value-of-an-integer-sequence
template<class T, T... Ints>
constexpr T get(std::integer_sequence<T, Ints...>, std::size_t i) {
    constexpr T arr[] = {Ints...};
    return arr[i];
}


template <std::size_t startNum, std::size_t... Length>
auto add_to_sequence(std::index_sequence<Length...>) {
  std::cout << "making preloop sequence. start num, number of numbers: " << startNum << "," << sizeof...(Length) << "\n";

  
  return std::index_sequence<(startNum + Length)...>{};




}


template <typename... TupleTemplate, std::size_t... indices>
auto tuple_slice(std::tuple<TupleTemplate...> knlTuple, std::index_sequence<indices...>) {

  auto newTuple = std::make_tuple(std::get<indices>(knlTuple)...);

  return newTuple;


}

template <typename... TupleTemplate, std::size_t... LoopNums>
auto chain_kernels(std::tuple<TupleTemplate...> knlTuple, Fuse<LoopNums...> fuseDirective) {
 
  
  std::cout << "Applying Fuse Transformation\n";

  std::cout << "Loop nums: ";
  ((std::cout << LoopNums << ", "), ...);
  std::cout << "\n";


  std::cout << "There are " << sizeof...(TupleTemplate) << " loops\n"; 
  std::cout << "Fusing " << sizeof...(LoopNums) << " of the loops\n";
  auto fuseNums = std::make_tuple(LoopNums...);

  constexpr auto firstNum = get(std::index_sequence<LoopNums...>(), 0);
  
  constexpr auto lastNum = get(std::index_sequence<LoopNums...>(), sizeof...(LoopNums) - 1);

  constexpr auto numPreLoop = firstNum;

  constexpr auto numPostLoop = sizeof...(TupleTemplate) - lastNum - 1;

  auto preFuseSequence = add_to_sequence<0>(std::make_index_sequence<numPreLoop>());
  auto postFuseSequence = add_to_sequence<lastNum+1>(std::make_index_sequence<numPostLoop>());
  std::cout << "firstNUm: " << firstNum << "\n";
  std::cout << "lastNum: " << lastNum << "\n";


  auto preFuseLoops = tuple_slice(knlTuple, preFuseSequence);
  auto postFuseLoops = tuple_slice(knlTuple, postFuseSequence);
  auto fuseLoops = std::make_tuple(std::get<LoopNums>(knlTuple)...);
  auto fusedLoops = shift_and_fuse_kernels(fuseLoops, std::make_index_sequence<sizeof...(LoopNums)>());

//  auto preFuseNums = 
  //auto postFuseNums 
  
  return 0;//return chain_kernels(knlTuple);
}

*/
/*
template <typename... TupleTemplate, std::size_t... LoopNums>
auto chain_kernels(std::tuple<TupleTemplate...> knlTuple, OverlappedTile<LoopNums...> fuseDirective) {
  std::cout << "Applying Overlapped Tiling Transformation\n";
 
  

   
  return chain_kernels(knlTuple);
} 
*/
/*
template <typename... TupleTemplate>
auto chain_kernels(std::tuple<TupleTemplate...> knlTuple) {
  std::cout << "Done applying transformations\n";

  return knlTuple;
} 


//MARK: Default chain codes


template <typename...Kernels>
void default_chain_kernels(Kernels&&... knls) {
  std::cout << "Performing default chain transformation pipeline\n";

  static auto kernelTuple = std::tuple(std::forward<Kernels>(knls)...);
  default_chain_shift_check_kernels(kernelTuple, std::index_sequence_for<Kernels...>{});

}//default_chain_kernels

template <typename...Kernels>
void default_chain_kernels() {
  std::cout << "Done with chain\n";
}
template<typename ...Kernels>
std::tuple<Kernels...> apply_shift_kernels(const std::vector<std::vector<int>> &, Kernels&&... knls) {
  std::cout << "TO IMPLEMENT: apply_shift_kernels\n";

  auto shiftedKnls = std::tuple(std::forward<Kernels>(knls)...);

  return shiftedKnls;
}

template<typename ...Kernels>
void execute_first_continue_kernels(auto knl1,  Kernels&&... knls) {

  knl1();
  default_chain_kernels(std::forward<Kernels>(knls)...);

}
template <typename ...Kernels>
int can_fuse_no_shift_kernels(Kernels&&...) {
  std::cout << "TO IMPLEMENT: can_fuse_no_shift_kernels\n";


  return 0; 
}

template <typename ...Kernels>
int can_fuse_with_shift_kernels(auto, Kernels&&...) {
  std::cout << "TO IMPLEMENT: can_fuse_with_shift_kernels\n";
  return 1;
}


template <typename Tuple, std::size_t... Is>
void default_chain_shift_check_kernels(const Tuple& knls, std::index_sequence<Is...>) {

  static int canFuseNoShift = can_fuse_no_shift_kernels(std::get<Is>(knls)...);
 
  if(canFuseNoShift) {
    std::cout << "Can fuse kernels without applying a shift\n";
  } else {
    std::cout << "Cannot fuse kernels without applying a shift\n";

    static std::vector<std::vector<int>> shiftAmounts = amount_to_shift_kernels(std::get<Is>(knls)...);
    static int canFuseWithShift = can_fuse_with_shift_kernels(shiftAmounts, std::get<Is>(knls)...);
    
    if(canFuseWithShift) {
      std::cout << "Can fuse kernels by applying a shift. Applying shift\n";
      auto shiftedKernelsTuple = apply_shift_kernels(shiftAmounts, std::get<Is>(knls)...);
      
      std::cout << "May want to execute pre-shift stuff here\n";
      default_chain_par_check_kernels(shiftedKernelsTuple, std::index_sequence<Is...>{});
      std::cout << "May want to execute post-shift shuff here\n";
      


    } else {
      std::cout << "Cannot fuse kernels by applying a shift. Executing first kernel and trying again\n";
      execute_first_continue_kernels(std::get<Is>(knls)...);
 
    }// canFuseWithShift
  }// canFuseNoShift




}//default chain_shift_check_kernels


template <typename ...Kernels>
int fusion_impedes_parallelism_kernels(Kernels&&...) {
  std::cout << "TO IMPLEMENT: fusion_impedes_parallelism_kernels\n";
  return 0;
}

template <typename ...Kernels>
int can_overlapped_tile_kernels(Kernels&&...) {
  std::cout << "TO IMPLEMENT: can_overlapped_tile_kernels\n";
  return 0;
}

template <typename ...Kernels>
int can_skew_kernels(Kernels&&...) {
  std::cout << "TO IMPLEMENT: can_skew_kernels\n";
  return 0;
}

template <typename ...Kernels>
void overlapped_tile_fuse_exec_kernels(Kernels&&...) {
  std::cout << "TO IMPLEMENT: overlapped_tile_fuse_exec_kernels\n";
  return;
}

template <typename ...Kernels>
void fuse_skew_exec_kernels(Kernels&&...) {
  std::cout << "TO IMPLEMENT: fuse_skew_exec_kernels\n";
  return;
}


template <typename ...Kernels>
void fuse_exec_kernels(Kernels&&...) {
  std::cout << "TO IMPLEMENT: fuse_exec_kernels\n";
  return;
}
template <typename Tuple, std::size_t... Is>
void default_chain_par_check_kernels(const Tuple& knls, std::index_sequence<Is...>) {

  std::cout << "Checking if fusion impedes parallelism\n";
  static int impedesParallelism = fusion_impedes_parallelism_kernels(std::get<Is>(knls)...);

  if(impedesParallelism) {
    std::cout << "Fusion impedes parallelism\n";

    static int canOverlappedTile = can_overlapped_tile_kernels(std::get<Is>(knls)...);

    if(canOverlappedTile) {
      std::cout << "Can apply overlapped tiling\n";
      overlapped_tile_fuse_exec_kernels(std::get<Is>(knls)...);
    } else {
      std::cout << "Cannot apply overlapped tiling\n";
      static int canSkew = can_skew_kernels(std::get<Is>(knls)...);
      if(canSkew) {
        std::cout << "Can apply skew\n";
        fuse_skew_exec_kernels(std::get<Is>(knls)...);
      } else {
        std::cout << "Cannot apply skew\n";
        execute_first_continue_kernels(std::get<Is>(knls)...);
      } // canSkew
    } //canOverlappedTile
  } else {
    std::cout << "Fusion does not impede parallelism\n";
   
    fuse_exec_kernels(std::get<Is>(knls)...);
  }//impedesParallelism


}//default_chain_par_check_kernels
*/
} //namespace RAJA
#endif
