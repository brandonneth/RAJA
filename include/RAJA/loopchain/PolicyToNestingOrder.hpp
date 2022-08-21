#ifndef RAJA_PolicyToNestingOrder_HPP
#define RAJA_PolicyToNestingOrder_HPP


namespace RAJA {


template <typename T>
struct policy_to_nesting_order_s {

};

template <typename InnerPolicy>
struct policy_to_nesting_order_s<KernelPolicy<InnerPolicy>> {

  auto operator()() {
    return policy_to_nesting_order_s<InnerPolicy>{}();
  }

};

template <camp::idx_t NestDepth, typename ExecPolicy, typename InnerPolicy>
struct policy_to_nesting_order_s<statement::For<NestDepth, ExecPolicy, InnerPolicy>> {
  std::vector<camp::idx_t> operator()() {
    std::vector<camp::idx_t> v = policy_to_nesting_order_s<InnerPolicy>{}();
    v.insert(v.begin(), NestDepth);

    return v;
  }
};

template <camp::idx_t LambdaId>
struct policy_to_nesting_order_s<statement::Lambda<LambdaId>> {
  auto operator()() {
    return std::vector<camp::idx_t>{};
  }
};

template <typename T> 
auto policy_to_nesting_order() {
    return policy_to_nesting_order_s<T>{}();
}
}
#endif