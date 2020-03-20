#ifndef RAJA_loopchain_kernel_HPP
#define RAJA_loopchain_kernel_HPP

#include "RAJA/config.hpp"

#include "RAJA/util/LoopChain.hpp"


namespace RAJA
{


 

template <typename PolicyType, typename SegmentTuple, typename... Bodies>
struct KernelW {

  const SegmentTuple &segments;
  std::tuple<Bodies...> bodies;
  const int numArgs;
  KernelW(SegmentTuple const & s, Bodies const &... b) : segments(s), bodies(b...), numArgs(camp::tuple_size<SegmentTuple>::value) {

  }

  KernelW(const KernelW & k) = default;
  KernelW(KernelW && k) = default;
  KernelW & operator=(const KernelW & k) = default;
  KernelW & operator=(KernelW && k) = default;


  std::vector<SymAccess> execute_symbolically() {

    std::cout << "Symbolically executing kernel, there are " << numArgs << " arguments to the kernel\n";
    SymIter i = SymIter("i0");
    SymIter j = SymIter("i1");

    auto func = std::get<0>(bodies);


    func(i,j);

    auto accesses1 = (i.accesses);
    auto accesses2 = (j.accesses);

    std::vector<SymAccess> allAccesses = std::vector<SymAccess>();

    for(long unsigned int i = 0; i < accesses1->size(); i++) {allAccesses.push_back(accesses1->at(i));}
    for(long unsigned int i = 0; i < accesses2->size(); i++) {allAccesses.push_back(accesses2->at(i));}
    return allAccesses;
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

  isl_union_map * deps = isl_union_map_union(raw, waw);
  deps = isl_union_map_union(deps, war);
  
  std::cout << "Dependences\n";
  p = isl_printer_print_union_map(p,deps);

  
  for(int i = 0; i < numKernels; i++) {
    for(int j = i + 1; j < numKernels; j++) {
      std::string srcLoop = "{L" + std::to_string(i) + "[" ;
      std::string dstLoop = "{L" + std::to_string(j) + "[";
      for(int k = 0; k < loopDims; k++) {
        srcLoop += "0";
        dstLoop += "n" + std::to_string(k);
        if( k != loopDims -1) {
          srcLoop += ",";
          dstLoop += ",";
        }
      }
      srcLoop += "] }";
      dstLoop += "] }";
      std::cout << "\nsource loop 0 iteration: " << srcLoop << "\n";
      std::cout << "dest loop iteration set: " << dstLoop << "\n";
      
      isl_union_set * srcSet = isl_union_set_read_from_str(ctx, srcLoop.c_str());
      isl_union_set * dstSet = isl_union_set_read_from_str(ctx, dstLoop.c_str());
      
      isl_union_set * fromSrc = isl_union_set_apply(srcSet, isl_union_map_copy(deps));
      
      isl_union_set * src2dst = isl_union_set_intersect(isl_union_set_copy(fromSrc), dstSet);

      std::cout << "\n deps from loop " << " i " << " to " << " j " << "\n";
      p = isl_printer_print_union_set(p, src2dst);
      
      // src2dst is a set of dependences from loop i to loop j. these dependences have loopDim dimensions
      // to convert these dependences to constraints, we have S_j_d - S_i_d = min(something, that i will figure out tomorrow.

    }
  }

  std::cout << "\n";
  return std::vector<std::vector<int>>();
}
  




} //namespace RAJA
#endif
