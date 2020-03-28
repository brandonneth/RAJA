#ifndef RAJA_loopchain_kernel_HPP
#define RAJA_loopchain_kernel_HPP

#include "RAJA/config.hpp"

#include "RAJA/util/LoopChain.hpp"


namespace RAJA
{


 

template<std::size_t s>
SymIter makeSymIter() {

  std::string iterName = "i" + std::to_string(s);
  
  return SymIter(iterName);

}

template <typename PolicyType, typename SegmentTuple, typename... Bodies>
struct KernelW {

  const SegmentTuple &segments;
  std::tuple<Bodies...> bodies;
  static constexpr int numArgs = camp::tuple_size<SegmentTuple>::value;
  KernelW(SegmentTuple const & s, Bodies const &... b) : segments(s), bodies(b...) {
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
    RAJA::kernel<PolicyType>(segments, std::get<Is>(bodies)...);

  }
  void operator () () {
     auto seq = std::index_sequence_for<Bodies...>{};
    execute(seq);
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

  std::cout << "Finished collecting access information, performing shift calculations\n";

  std::cout << "Gathered access information for " << accessLists.size() << " kernels, each with " << numArgs << " dims\n";

  auto numKernels = accessLists.size();

  auto loopDims = numArgs;


  auto readString = read_string(accessLists, loopDims);
  auto writeString = write_string(accessLists, loopDims); 


  std::cout << "read string: "<< readString << "\n";
  std::cout << "write string: " << writeString << "\n";

  isl_ctx* ctx = isl_ctx_alloc();
  isl_printer * p = isl_printer_to_file(ctx, stdout);

  isl_union_map * reads = isl_union_map_read_from_str(ctx, readString.c_str());
  isl_union_map * writes = isl_union_map_read_from_str(ctx, writeString.c_str());

  isl_union_map * reads_inverse = isl_union_map_reverse(isl_union_map_copy(reads));
  isl_union_map * writes_inverse = isl_union_map_reverse(isl_union_map_copy(writes));

  isl_union_map * raw = isl_union_map_apply_range(isl_union_map_copy(writes), isl_union_map_copy(reads_inverse));
  isl_union_map * waw = isl_union_map_apply_range(isl_union_map_copy(writes), isl_union_map_copy(writes_inverse));
  isl_union_map * war = isl_union_map_apply_range(isl_union_map_copy(reads), isl_union_map_copy(writes_inverse));

  
  std::cout << "\nwaw\n";
  p = isl_printer_print_union_map(p,waw);
  std::cout << "\nraw\n";
  p = isl_printer_print_union_map(p,raw);
  std::cout << "\nwar\n";
  p = isl_printer_print_union_map(p,war);



  isl_union_map * deps = isl_union_map_union(raw, waw);
  deps = isl_union_map_union(deps, war);
  
  std::cout << "\nDependences\n";
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
  std::cout << "\nshift set" << shiftSet << "\n";
 
  std::string allPossibleShiftsString = "{" + shiftSet + ":" + nonNegConstraint + "}";
  std::cout << "\nnonneg shift set string: " << allPossibleShiftsString << "\n";
  isl_union_set * allPossibleShifts = isl_union_set_read_from_str(ctx, allPossibleShiftsString.c_str());

  std::cout << "\nallPossibleShifts\n";

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
      std::cout << "\nsource loop: " << srcLoop << "\n";
      std::cout << "dest loop: " << dstLoop << "\n";
    
      isl_union_set * srcSet = isl_union_set_read_from_str(ctx, srcLoop.c_str());
      isl_union_set * dstSet = isl_union_set_read_from_str(ctx, dstLoop.c_str());
      std::cout << "srcSet\n";
      p= isl_printer_print_union_set(p,srcSet);

      std::cout << "\ndstSet\n";
      p = isl_printer_print_union_set(p,dstSet);
      isl_union_set * depsFromSrc = isl_union_set_apply(srcSet, isl_union_map_copy(deps));
      
      std::cout << "\ndepsFromSrc\n";
      p = isl_printer_print_union_set(p,depsFromSrc);


      isl_union_set * src2dstUnionSet = isl_union_set_intersect(isl_union_set_copy(depsFromSrc), dstSet);
      std::cout << "\nsrc2dstUnionSet\n";
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
        std::cout << "shift constraint: " << constraint << "\n";
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
  



template <int tileDimension, typename ExecPol, typename Container, typename... Bodies>
auto tile_kernel(KernelW<KernelPolicy<ExecPol>,Container,Bodies...> & knl, int numTiles) {
  
  std::cout << "tiling kernel with " << knl.numArgs <<"on dimension" << tileDimension <<  "\n"; 
  using TiledPol = RAJA::KernelPolicy<
    RAJA::statement::For<knl.numArgs, RAJA::seq_exec,ExecPol>
  >;
  
  auto indices = camp::make_idx_seq_t<knl.numArgs>();
  auto newLambda = [=] (auto i, auto j, auto iTile) {
    std::cout << "tiling lambda call: "<< i << " "<< j << " " << iTile << "\n";
  };

  //auto newTuple = RAJA::make_tuple(iSegment,jSegment,tileSegment);

  return knl;
  //return makeKernel<TiledPol>(newTuple, newLambda);
}


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
 


  
  auto tiledSegment = RAJA::RangeSegment(0,tileSize);

  using OVERLAPPEDPOL = RAJA::KernelPolicy<
    RAJA::statement::For<0, RAJA::seq_exec, //tile loop
      RAJA::statement::For<Is+1, RAJA::seq_exec, //overlap1 loop
        RAJA::statement::Lambda<Is, RAJA::statement::Segs<Is+1,0>>>...,
      RAJA::statement::For<4, RAJA::seq_exec, // tiled loop part
        RAJA::statement::Lambda<Is+sizeof...(Kernels), RAJA::statement::Segs<4,0>>...
      >
    >
  >;

  auto totalTuple = RAJA::make_tuple(tileSegment,std::get<Is>(overlappedSegments)..., tiledSegment);

  auto tiledKnl = RAJA::makeKernel<OVERLAPPEDPOL>(totalTuple, std::get<Is>(overlapLambdas)...,std::get<Is>(overlapLambdas)...);
 
  return tiledKnl; 
}





template <typename ...Kernels> 
auto overlapped_tiling_lambda_nofirst(Kernels&&... knls) {

  int numTiles = 4;
  
  int upperBound = 16;


  auto otl = [=] (auto t, auto i) {
    auto new_i = t * upperBound / numTiles + i;
  };



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

} //namespace RAJA
#endif
