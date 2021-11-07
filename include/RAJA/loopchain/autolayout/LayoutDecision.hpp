#ifndef RAJA_LayoutDecision_HPP
#define RAJA_LayoutDecision_HPP

#include "RAJA/config.hpp"


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


  template <camp::idx_t I, typename ViewType, typename LayoutType, typename MomentView>
  int check_and_update_layout(LayoutTypes & moment, MomentView & moment_view, ViewType & view_to_change, LayoutType new_layout) {
    if constexpr (std::is_same<ViewType, MomentView>::value) {
     
      if(&moment_view == &view_to_change) {
        std::cout << "The views are refreences to the same views\n";
        get<I>(moment) = new_layout;
      }
      return 0;
    } else {
      return 0;
    }
  }

  template <typename ViewType, typename LayoutType, camp::idx_t...Is>
  void update_moment(LayoutTypes & moment, ViewType & view, LayoutType layout, camp::idx_seq<Is...>) {
    auto results = tuple(check_and_update_layout<Is>(moment, get<Is>(views), view, layout)...);
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

}; //LayoutDecisions

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
