

#ifndef RAJA_LoopChain_HPP
#define RAJA_LoopChain_HPP

#include "RAJA/config.hpp"

#include "RAJA/loopchain/Utils.hpp"
#include "RAJA/loopchain/SymExec.hpp"
#include "RAJA/loopchain/KernelWrapper.hpp"
#include "RAJA/loopchain/Chain.hpp"



#include "RAJA/loopchain/Transformations.hpp"
#include "RAJA/loopchain/ISLAnalysis.hpp"
#include "RAJA/loopchain/KernelConversion.hpp"

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

#endif
