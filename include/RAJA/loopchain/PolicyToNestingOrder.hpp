#ifndef RAJA_PolicyToNestingOrder_HPP
#define RAJA_PolicyToNestingOrder_HPP
#include "RAJA/config.hpp"
#include "RAJA/util/resource.hpp"
#include "RAJA/util/basic_mempool.hpp"
#include "RAJA/pattern/forall.hpp"
#include "RAJA/pattern/kernel.hpp"
namespace RAJA {


template <typename T, camp::idx_t SearchIdx>
struct policy_to_nesting_order_s {

};

template <typename InnerPolicy, camp::idx_t SearchIdx>
struct policy_to_nesting_order_s<RAJA::KernelPolicy<InnerPolicy>, SearchIdx> {

  auto operator()() {
    return policy_to_nesting_order_s<InnerPolicy, SearchIdx>{}();
  }

};

template <camp::idx_t NestDepth, typename ExecPolicy, typename...InnerPolicy, camp::idx_t SearchIdx>
struct policy_to_nesting_order_s<statement::For<NestDepth, ExecPolicy, InnerPolicy...>, SearchIdx> {
  std::vector<camp::idx_t> operator()() {
    std::vector<camp::idx_t> v = policy_to_nesting_order_s<camp::tuple<InnerPolicy...>, SearchIdx>{}();
    v.insert(v.begin(), NestDepth);

    return v;
  }
};

template <camp::idx_t LambdaId, typename...LambdaOther, camp::idx_t SearchIdx>
struct policy_to_nesting_order_s<statement::Lambda<LambdaId, LambdaOther...>, SearchIdx> {
  auto operator()() {
    if constexpr (LambdaId == SearchIdx) {
      return std::vector<camp::idx_t>{};
    } else {
      return std::vector<camp::idx_t>{{-1}};
    }
  }
};

template <typename Statement1, typename...Statements, camp::idx_t SearchIdx>
struct policy_to_nesting_order_s<camp::tuple<Statement1, Statements...>, SearchIdx> {
  auto operator()() {
    std::vector<camp::idx_t> statement1Piece = policy_to_nesting_order_s<Statement1, SearchIdx>{}();
    if (statement1Piece.size() == 0 || statement1Piece.back() != -1) {
      return statement1Piece;
    } else {
      if constexpr (sizeof...(Statements) != 0) {
        return policy_to_nesting_order_s<camp::tuple<Statements...>, SearchIdx>{}();
      } else {
        return std::vector<camp::idx_t>{{-1}};
      }
    }
  }
};

template <typename T, idx_t SearchIdx> 
auto policy_to_nesting_order() {
    return policy_to_nesting_order_s<T, SearchIdx>{}();
}

template <typename T> 
auto policy_to_nesting_order() {
    return policy_to_nesting_order<T,0>();
}

}
#endif
