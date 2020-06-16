

#ifndef RAJA_KernelWrapper_HPP
#define RAJA_KernelWrapper_HPP

#include "RAJA/config.hpp"
#include "RAJA/loopchain/SymExec.hpp"
#include "RAJA/pattern/kernel.hpp"

#include <vector>
#include <string>
namespace RAJA
{



template <typename KernelPol, typename SegmentTuple, typename... Bodies>
struct KernelWrapper {

  const SegmentTuple segments;
  const camp::tuple<Bodies...> bodies;
 
  std::vector<camp::idx_t> overlapAmounts;
  std::vector<camp::idx_t> tileSizes;
  static constexpr int numArgs = camp::tuple_size<SegmentTuple>::value;

  KernelWrapper(SegmentTuple  _segments, const Bodies&... _bodies) : 
    segments(_segments), bodies(_bodies...) {
     overlapAmounts = std::vector<camp::idx_t>();
     tileSizes = std::vector<camp::idx_t>();
  }
  
  KernelWrapper(const KernelWrapper &) = default;
  KernelWrapper(KernelWrapper &&) = default;
  KernelWrapper & operator=(const KernelWrapper &) = default;
  KernelWrapper & operator=(KernelWrapper &&) = default;


  template <std::size_t Idx>
  SymIterator make_sym_iterator() {
    std::string iteratorName = "i" + std::to_string(Idx);
    return SymIterator(iteratorName);
  }

  template <std::size_t... Is>
  auto make_iterator_tuple(camp::idx_seq<Is...>) {
    auto iterators = camp::make_tuple((make_sym_iterator<Is>())...);
    return iterators;
  }
  
  std::vector<SymAccess> collect_accesses() {
    return std::vector<SymAccess>();
  }

  template <typename... Iterators>
  std::vector<SymAccess> collect_accesses(SymIterator iterator, Iterators&&... rest) {
    std::vector<SymAccess> accesses = collect_accesses(std::forward<Iterators>(rest)...);

    for(long unsigned int i = 0; i < iterator.accesses->size(); i++) {
      accesses.push_back(iterator.accesses->at(i));
    }
 
    return accesses;
  }

  template <std::size_t... Is>
  auto collect_accesses_from_iterators(auto iterators, camp::idx_seq<Is...>) {
    return collect_accesses(camp::get<Is>(iterators)...);
  }
  
  std::vector<SymAccess> execute_symbolically() {
    auto iterators = make_iterator_tuple(camp::make_idx_seq_t<numArgs>());

    auto func = camp::get<0>(bodies);

    camp::invoke(iterators, func);

    auto accesses = collect_accesses_from_iterators(iterators, camp::make_idx_seq_t<numArgs>());
 
    return accesses;
  }

  template <camp::idx_t... Is>
  void execute(camp::idx_seq<Is...>) const { 
  
    util::PluginContext context{util::make_context<KernelPol>()};
    util::callPreLaunchPlugins(context);

    using segment_tuple_t = typename IterableWrapperTuple<camp::decay<SegmentTuple>>::type;

    auto params = RAJA::make_tuple();
    using param_tuple_t = camp::decay<decltype(params)>;

    if(overlapAmounts.size() != 0 && tileSizes.size() != 0) {
      using loop_data_t = internal::LoopData<KernelPol, segment_tuple_t, param_tuple_t, camp::decay<Bodies>...>;

      loop_data_t loop_data(overlapAmounts, tileSizes, make_wrapped_tuple(segments), params, camp::get<Is>(bodies)...);

      RAJA_FORCEINLINE_RECURSIVE
      internal::execute_statement_list<KernelPol>(loop_data);

      util::callPostLaunchPlugins(context);
    } else {
      RAJA::kernel<KernelPol>(segments, camp::get<Is>(bodies)...);
    }







  } //execute

  void operator() () const {
    auto seq = camp::make_idx_seq_t<sizeof...(Bodies)>{};
    execute(seq);
  }
}; // KernelWrapper

template <typename KernelPol, typename SegmentTuple, typename... Bodies>
KernelWrapper<KernelPol, SegmentTuple, Bodies...> 
make_kernel(const SegmentTuple & segment,   Bodies const &... bodies) {

  return KernelWrapper<KernelPol,SegmentTuple,Bodies...>(segment, bodies...);
}


template <typename ExecPol, typename Segment, typename Body> 
auto make_forall(Segment segment, const Body & body) {
  using KernPol = 
    RAJA::KernelPolicy<
      statement::For<0,ExecPol,
        statement::Lambda<0>
      >
    >;     
  
  return KernelWrapper<KernPol, camp::tuple<Segment>, Body>(camp::make_tuple(segment), body);
}


} //namespace RAJA


#endif

