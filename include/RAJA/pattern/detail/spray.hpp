#ifndef RAJA_PATTERN_DETAIL_REDUCE_SPRAY_HPP
#define RAJA_PATTERN_DETAIL_REDUCE_SPRAY_HPP

namespace RAJA {


template<typename ReductionPolicy, typename ContentType>
class SprayReduction {
public:
  mutable ContentType * content;

  virtual ContentType& operator() (long i) const
  {
    std::cout << "Base class accessor\n";
    return content[i];
  }



};
}
#include "RAJA/pattern/detail/spray/dense.hpp"
#include "RAJA/pattern/detail/spray/atomic.hpp"
#include "RAJA/pattern/detail/spray/map.hpp"
#include "RAJA/pattern/detail/spray/block.hpp"
#include "RAJA/pattern/detail/spray/keeper.hpp"
#endif
