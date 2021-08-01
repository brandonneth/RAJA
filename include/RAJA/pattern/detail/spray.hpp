#ifndef RAJA_PATTERN_DETAIL_REDUCE_SPRAY_HPP
#define RAJA_PATTERN_DETAIL_REDUCE_SPRAY_HPP

template <typename ReductionPolicy, typename ContentType>
class SprayReduction
{
  virtual ContentType &operator()(int idx) = 0;
}


#include "RAJA/pattern/detail/spray/dense.hpp"
#include "RAJA/pattern/detail/spray/atomic.hpp"
#include "RAJA/pattern/detail/spray/map.hpp"
#include "RAJA/pattern/detail/spray/block.hpp"
#include "RAJA/pattern/detail/spray/keeper.hpp"
#endif
