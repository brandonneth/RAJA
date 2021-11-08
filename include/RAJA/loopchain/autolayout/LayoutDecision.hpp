#ifndef RAJA_LayoutDecision_HPP
#define RAJA_LayoutDecision_HPP

#include "RAJA/config.hpp"
#include "RAJA/loopchain/Utils.hpp"

namespace RAJA {

template <typename ViewType, typename LayoutType>
auto permute_view(ViewType & v, LayoutType l);


template <camp::idx_t NumKernels, typename...Views>
struct LayoutDecisions {

  constexpr static camp::idx_t NumViews = sizeof...(Views);

  using LayoutTypes = camp::tuple<typename std::remove_reference<Views>::type::layout_type...>;
   
  std::array<LayoutTypes, NumKernels + 1> moments;

  camp::tuple<Views...> views;
  LayoutDecisions(camp::tuple<Views...> _views) :views(_views) {
    
    std::cout << "instatiating a LayoutDecisions object with " << NumViews << " views\n";
   
  }

  LayoutDecisions(const LayoutDecisions & other) : moments(other.moments), views(other.views) {
    
  }
/*
  LayoutDecisions& operator=(LayoutDecisions other) {
    moments = other.moments;
    views = other.views;
  }*/

  template <camp::idx_t I, typename ViewType, typename LayoutType, typename MomentView>
  int check_and_update_layout(LayoutTypes & moment, MomentView & moment_view, ViewType & view_to_change, LayoutType new_layout) {
    if constexpr (std::is_same<ViewType, MomentView>::value) {
     
      if(&moment_view == &view_to_change) {
        std::cout << "The views are refreences to the same views\n";
        get<I>(moment) = new_layout;
        return 1;
      }
      return 0;
    } else {
      return 0;
    }
  }

  template <typename ViewType, typename LayoutType, camp::idx_t...Is>
  void update_moment(LayoutTypes & moment, ViewType & view, LayoutType layout, camp::idx_seq<Is...>) {
    auto results = (check_and_update_layout<Is>(moment, get<Is>(views), view, layout) + ...);
    if(results == 0) {
      std::cout << "call to update_moment did not change any of the layouts\n";
    }
  }

  template <typename ViewType, typename LayoutType>
  void add_decision(ViewType & view_to_change, LayoutType new_layout, camp::idx_t when) {
    update_moment(moments[when], view_to_change, new_layout, idx_seq_for(moments[when]));
  }




  template <camp::idx_t...Is>
  auto kernels_helper(LayoutTypes moment, camp::idx_seq<Is...>) {
    return make_tuple(permute_view(get<Is>(views), get<Is>(moment))...);
  }

  auto kernels(camp::idx_t i) {
    auto moment = moments[i];

    auto knls = kernels_helper(moment, idx_seq_for(moment));
    return knls;
  }

  template <camp::idx_t NumKernels_, typename...Views_>
  friend std::ostream & operator<<(std::ostream& os, const LayoutDecisions & ld);
}; //LayoutDecisions

template <typename Moment, camp::idx_t I, camp::idx_t...Is>
void output_moment(std::ostream&os, const Moment & moment, camp::idx_seq<I, Is...>) {
  auto layout = get<I>(moment);
  os << "_view" << I << "_" ;
  auto perm = layout_to_perm(layout);

  for(int i = 0; i < perm.size(); i++) {
  os << perm[i] << "_";
  }
  os << "_";
  if constexpr (sizeof...(Is) == 0) { 
    return;
  } else {
    output_moment(os, moment, camp::idx_seq<Is...>{});
  }
}


template <camp::idx_t CurrKernel, camp::idx_t NumKernels, typename...Views, camp::idx_t...Is>
void ld_output_helper(std::ostream&os, const LayoutDecisions<NumKernels,Views...> & ld, camp::idx_seq<Is...> seq) {
  auto moment = ld.moments[CurrKernel];
  output_moment(os, moment, seq);
  if constexpr (CurrKernel == NumKernels) {
    return;
  } else {
    os <<"_" << "knl" << CurrKernel << "_";
    ld_output_helper<CurrKernel+1, NumKernels>(os, ld, seq);
    return;
  }
}

template <camp::idx_t NumKernels, typename...Views>
std::ostream & operator<<(std::ostream& os, const LayoutDecisions<NumKernels,Views...> & ld) {
  ld_output_helper<0,NumKernels>(os, ld, idx_seq_from_to<0,sizeof...(Views)>());
  return os;
} 


template <typename...Views, typename KernelType1, typename KernelType2, camp::idx_t...Is>
auto decide_moment(int moment_number, LayoutDecisions<2,Views...> & decisions,  KernelType1 knl1, KernelType2 knl2, camp::idx_seq<Is...>) {

  
}

template <typename...Views, typename KernelType1, typename KernelType2>
auto make_decisions(LayoutDecisions<2,Views...> & decisions,  KernelType1 knl1, KernelType2 knl2) {
  for(int i = 0; i < 3; i++) {
    decide_moment(i, decisions, knl1, knl2, idx_seq_from_to<0,sizeof...(Views)>());
  }
}
template <typename...Views, typename KernelType1, typename KernelType2>
auto fill_remaining_decisions(LayoutDecisions<2,Views...> & decisions,  KernelType1 knl1, KernelType2 knl2) {

  return;
}

template <typename...Views, typename KernelType1, typename KernelType2>
auto manipulate_layouts(LayoutDecisions<2,Views...> & decisions,  KernelType1 knl1, KernelType2 knl2) {

  fill_remaining_decisions(decisions, knl1, knl2);

  auto conv0 = decisions.kernels(0);
  auto conv1 = decisions.kernels(1);
  auto conv2 = decisions.kernels(2);

  auto all_knls = tuple_cat(conv0, make_tuple(knl1), conv1, make_tuple(knl2), conv2);

  return chain(all_knls);
}


template <camp::idx_t NumKernels, typename... ViewTypes>
auto create_layout_decisions(tuple<ViewTypes...> views) {
  using DecisionType = LayoutDecisions<NumKernels, ViewTypes... >;

  return DecisionType(views);
}


} //namespace RAJA 


#endif
