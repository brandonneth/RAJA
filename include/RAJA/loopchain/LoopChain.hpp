

#ifndef RAJA_LoopChain_HPP
#define RAJA_LoopChain_HPP

#include "RAJA/config.hpp"

#include "RAJA/loopchain/Utils.hpp"
#include "RAJA/loopchain/SymExec.hpp"
#include "RAJA/loopchain/KernelWrapper.hpp"
#include "RAJA/loopchain/Transformations.hpp"
#include "RAJA/loopchain/ISLAnalysis.hpp"
#include "RAJA/loopchain/KernelConversion.hpp"


//
// API Definitions
//
namespace RAJA
{

template <typename...Ts>
auto shift(auto knl, tuple<Ts...> shiftAmountTuple);

template <typename...Ts>
auto shift(auto knl, Ts... shiftAmounts);

template <typename FusedKernelType, typename...Knls>
FusedKernelType fuse(Knls...knls);

template <typename FusedKernelType, typename...Knls>
FusedKernelType shift_and_fuse(Knls...knls);

template <typename...Knls>
auto overlapped_tile_no_fuse(Knls...knls);

template <typename...Knls>
auto overlapped_tile_fuse(Knls...knls);

template <typename...Knls>
auto chain(Knls...knls);

} //namespace RAJA


//
// Implementation
//
namespace RAJA 
{

template <camp::idx_t mainKnlIdx, typename... KernelTupleType>
struct FusedKernels {

  const camp::tuple<KernelTupleType...> knlTuple;

  static constexpr  camp::idx_t numKnls = sizeof...(KernelTupleType);

  FusedKernels(const auto & knlTuple_) : knlTuple(knlTuple_) {};
  
  auto pre_knls() {
    return slice_tuple<0,mainKnlIdx>(knlTuple);
  }

  auto main_knl() {
    return camp::get<mainKnlIdx>(knlTuple);
  }

  auto post_knls() {
    return slice_tuple<mainKnlIdx+1, numKnls>(knlTuple);
  }


  template <camp::idx_t I>
  RAJA_INLINE
  void execute(camp::idx_seq<I>) {
    camp::get<I>(knlTuple)();
  }
  template <camp::idx_t I, camp::idx_t... Is>
  RAJA_INLINE
  void execute(camp::idx_seq<I, Is...>) {
    camp::get<I>(knlTuple)();
    execute(camp::idx_seq<Is...>{});
  }
  
  RAJA_INLINE
  void operator() () {
    auto seq = camp::make_idx_seq_t<sizeof...(KernelTupleType)>{};
    execute(seq); 
  }

};

template <camp::idx_t MainKnlIdx, typename...KnlTupleType>
auto fused_kernels(camp::tuple<KnlTupleType...> knlTuple) {
  return FusedKernels<MainKnlIdx,KnlTupleType...>(knlTuple);
} //fused_kernels


template <typename... KernelTupleType>
struct LoopChain {

  const camp::tuple<KernelTupleType...> knlTuple;
   
  LoopChain(const auto & knlTuple_) : knlTuple(knlTuple_) {}
  
  
  template <camp::idx_t I>
  void execute(camp::idx_seq<I>) {
    camp::get<I>(knlTuple)();
  }
  template <camp::idx_t I, camp::idx_t... Is>
  void execute(camp::idx_seq<I, Is...>) {
    camp::get<I>(knlTuple)();
    execute(camp::idx_seq<Is...>{});
  }
  
  void operator() () {
    auto seq = camp::make_idx_seq_t<sizeof...(KernelTupleType)>{};
    execute(seq); 
  }

}; //LoopChain


template <typename...Knls>
auto loop_chain(camp::tuple<Knls...> knlTuple) {
  return LoopChain<Knls...>(knlTuple);
}
template <typename...Knls>
auto loop_chain(Knls...knls) {
  auto knlTuple = make_tuple(knls...);
  return LoopChain<Knls...>(knlTuple);
} //loop_chain


/*

template <typename... KernelTupleTemplate, 
          typename... TransformationTupleTemplate>
auto chain(camp::tuple<KernelTupleTemplate...> knlTuple,
           camp::tuple<TransformationTupleTemplate...> transformationTuple) {
  return LoopChain<KernelTupleTemplate...>(knlTuple);
}

template <typename... KernelTupleTemplate, 
          typename... TransformationTupleTemplate,
          typename... KernelTemplate, 
          typename... Rest>
auto chain(const camp::tuple<KernelTupleTemplate...> & knlTuple, 
           camp::tuple<TransformationTupleTemplate...> transformationTuple,
           const KernelWrapper<KernelTemplate...> & knl, 
           Rest&&... rest) {
  
  auto newTuple = camp::tuple_cat_pair(knlTuple, camp::make_tuple(knl));
  
  return chain(newTuple, transformationTuple, std::forward<Rest>(rest)...);

}

template <typename...KernelTupleTemplate, 
          typename...TransformationTupleTemplate, 
          camp::idx_t LoopId,
          typename... ShiftAmounts,
          typename... Rest>
auto chain(const camp::tuple<KernelTupleTemplate...> & knlTuple,
           camp::tuple<TransformationTupleTemplate...> transformationTuple,
           Shift<LoopId,ShiftAmounts...> shift,
           Rest&&... rest) {

  auto knlToShift = camp::get<LoopId>(knlTuple);
  auto shiftedKnl = shift_kernel(knlToShift, shift);

  auto newKnlTuple = replace_index_with<LoopId>(knlTuple, make_tuple(shiftedKnl));
  
  return chain(newKnlTuple, transformationTuple, std::forward<Rest>(rest)...);
}
 
template <typename...KernelTupleTemplate,
          typename...TransformationTupleTemplate,
          camp::idx_t LoopId1, camp::idx_t LoopId2,
          typename... Rest>
auto chain(const camp::tuple<KernelTupleTemplate...> & knlTuple,
           camp::tuple<TransformationTupleTemplate...> transformationTuple,
           Fusion<LoopId1, LoopId2> fuse,
           Rest&&... rest) {


  auto fusedKnl = fuse_kernels(camp::get<LoopId1>(knlTuple), camp::get<LoopId2>(knlTuple));

   
  auto tmpKnlTuple = remove_index<LoopId2>(knlTuple);
  auto newKnlTuple = replace_index_with<LoopId1>(tmpKnlTuple, fusedKnl);

  return chain(newKnlTuple, transformationTuple, std::forward<Rest>(rest)...);
}

template <typename...KernelTupleTemplate,
          typename...TransformationTupleTemplate,
          camp::idx_t LoopId, typename...OverlappedTileType,
          typename... Rest>
auto chain(const camp::tuple<KernelTupleTemplate...> & knlTuple,
           camp::tuple<TransformationTupleTemplate...> transformationTuple,
           OverlappedTile<LoopId, OverlappedTileType...> otile,
           Rest&&... rest) {


  auto tiledKnl = overlapped_tile_kernel(camp::get<LoopId>(knlTuple), otile.overlapAmounts, otile.tileSizes);

   
  auto newKnlTuple = replace_index_with<LoopId>(knlTuple, make_tuple(tiledKnl));

  return chain(newKnlTuple, transformationTuple, std::forward<Rest>(rest)...);
}


template <typename... KernelTemplate, typename... Rest>
auto chain(const KernelWrapper<KernelTemplate...> & knl, Rest&&...  rest) {

  auto knlTuple = RAJA::make_tuple(knl);

  return chain(knlTuple, RAJA::make_tuple(), std::forward<Rest>(rest)...);
}

*/

} //namespace RAJA

#endif
