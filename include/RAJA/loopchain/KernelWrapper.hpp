// Contains code for the kernel wrapper objects and the functions that create them.

#ifndef RAJA_KernelWrapper_HPP
#define RAJA_KernelWrapper_HPP

#include "RAJA/config.hpp"
#include "RAJA/loopchain/SymExec.hpp"


#include <vector>
#include <string>
namespace RAJA
{



template <typename ExecPol, typename SegmentTuple, typename... Bodies>
struct KernelWrapper {

  const SegmentTuple &segments;
  camp::tuple<Bodies...> bodies;
  
  static constexpr int numArgs = camp::tuple_size<SegmentTuple>::value;

  KernelWrapper(SegmentTuple const & _segments, Bodies const &... _bodies) : 
    segments(_segments), bodies(_bodies...) {
     
  }
  
  KernelWrapper(const KernelWrapper &) = default;
  KernelWrapper(KernelWrapper &&) = default;
  KernelW & operator=(const KernelWrapper &) = default;
  KernelWrapper & operator=(KernelWrapper &&) = default;


  template <std::size_t Idx>
  SymIterator make_sym_iterator() {
    std::string iteratorName = "i" + std::to_string(Idx);
    return SymIterator(iteratorName);
  }

  template <std::size_t... Is>
  auto make_iterator_tuple(camp::index_sequence<Is...>) {
    auto iterators = camp::tuple((make_sym_iterator<Is>())...);
    return iterators;
  }
  
  std::vector<SymAccess> collect_accesses() {
    return std::vector<SymAccess>();
  }

  template <typename... Iterators>
  std::vector<SymAccess> collect_accesses(SymIterator iterator, Iterators&&... rest) {
    std::vector<SymAccess accesses = collect_accesses(std::forward<Iterators>(rest)...);

    for(long unsigned int i = 0; i < iterator.accesses->size(); i++) {
      accesses.push_back(iterator.accesses->at(i));
    }
 
    return allAccesses;
  }

  template <std::size_t... Is>
  auto collect_accesses_from_iterators(auto iterators, camp::index_sequence<Is...>) {
    return collect_accesses(camp::get<Is>(iterators)...);
  }
  
  std::vector<SymAccess> execute_symbolically() {
    auto iterators = make_iterator_tuple(camp::make_index_sequence<numArgs>());

    auto func = camp::get<0>(bodies);

    camp::apply(func, iterators);

    auto accesses = collect_accesses_from_iterators(iterators, camp::make_index_sequence<numArgs>());
  }

  template <std::size_t... Is>
  void execute(camp::index_sequence<Is...>) {
  
    util::PluginContext context{util::make_context<PolicyType>()};
    util::callPreLaunchPlugins(context);

    using segment_tuple_t = typename IterableWrapperTuple<camp::decay<SegmentTuple>>::type;

    auto params = RAJA::make_tuple();
    using param_tuple_t = camp::decay<decltype(params)>;

    //For when i introduce overlapped tiling
    /*
    if(overlaps.size() != 0){
      using loop_data_t = internal::LoopData<PolicyType, segment_tuple_t, param_tuple_t, camp::decay<Bodies>...>;
      loop_data_t loop_data(overlaps, make_wrapped_tuple(segments), params, std::get<Is>(bodies)...);

      RAJA_FORCEINLINE_RECURSIVE
      internal::execute_statement_list<PolicyType>(loop_data);

    } else {
      RAJA::kernel<PolicyType>(segments, std::get<Is>(bodies)...);
    }
    */

    RAJA::kernel<ExecPol>(segments, camp::get<Is>(bodies)...);
  }

  void operator() () {
    auto seq = camp::index_sequence_for<Bodies...>{};
    execute(seq);
  }
}


} //namespace RAJA


#endif

