#ifndef RAJA_AutoLayout_HPP
#define RAJA_AutoLayout_HPP

#include "RAJA/config.hpp"

namespace RAJA {

template <typename T>
struct policy_to_nesting_order {

};

template <typename InnerPolicy>
struct policy_to_nesting_order<KernelPolicy<InnerPolicy>> {

  auto operator()() {
    return policy_to_nesting_order<InnerPolicy>{}();
  }

};

template <camp::idx_t NestDepth, typename ExecPolicy, typename InnerPolicy>
struct policy_to_nesting_order<statement::For<NestDepth, ExecPolicy, InnerPolicy>> {
  std::vector<camp::idx_t> operator()() {
    std::vector<camp::idx_t> v = policy_to_nesting_order<InnerPolicy>{}();
    v.insert(v.begin(), NestDepth);

    return v;
  }
};

template <camp::idx_t LambdaId>
struct policy_to_nesting_order<statement::Lambda<LambdaId>> {
  auto operator()() {
    return std::vector<camp::idx_t>{};
  }
};

}


#endif
