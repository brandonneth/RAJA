#ifndef RAJA_AutoLayout_HPP
#define RAJA_AutoLayout_HPP

#include "RAJA/config.hpp"

namespace RAJA {

template <typename T>
struct policy_to_argument_order {

};

template <typename InnerPolicy>
struct policy_to_argument_order<KernelPolicy<InnerPolicy>> {

  auto operator()() {
    return std::vector<camp::idx_t>();
  }

};



}


#endif
