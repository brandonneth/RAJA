#ifndef RAJA_LayoutDecision_HPP
#define RAJA_LayoutDecision_HPP

#include "RAJA/config.hpp"
#include "RAJA/loopchain/Utils.hpp"

namespace RAJA {

template <typename ViewType, typename LayoutType>
auto permute_view(ViewType & v, LayoutType l);



template <camp::idx_t NumKernels, typename...ViewRefs>
struct FormatDecisions {
  constexpr static camp::idx_t NumViews = sizeof...(ViewRefs);
  constexpr static camp::idx_t NumMoments = NumKernels+1;
  using LayoutTupleType = camp::tuple<typename std::remove_reference<ViewRefs>::type::layout_type...>;

  
  camp::tuple<ViewRefs...> views;
  std::array<std::array<bool, NumViews>, NumMoments> user_decided;

  std::array<LayoutTupleType, NumMoments> moments;

  LayoutTupleType input_formats;

  template <camp::idx_t...Is>
  FormatDecisions(camp::tuple<ViewRefs...> _views, camp::idx_seq<Is...>) : views(_views) {
    input_formats = make_tuple(get<Is>(views).get_layout()...);

    for(int i0 = 0; i0 < NumMoments; i0++ ) {
      for(auto i1 = 0; i1 < NumViews; i1++) {
        user_decided[i0][i1] = false;
      }
    }
  }

  

  template <camp::idx_t I, typename ViewType, typename InternalViewType>
  camp::idx_t view_index_helper(ViewType & view, InternalViewType & internal_view) {
   
    if constexpr (std::is_same<ViewType, InternalViewType>::value) {
      if (&view == &internal_view) {
        return I;
      } else {
        return -1;
      }
    } else {
      return -1;
    }
  }

  template <typename ViewType, camp::idx_t...Is>
  camp::idx_t view_index(ViewType & view, camp::idx_seq<Is...>) {
    camp::idx_t index = vmax(view_index_helper<Is>(view, get<Is>(views))...);
    return index;
  }

  template <typename ViewType, typename LayoutType>
  bool change_to_before(ViewType & view, LayoutType layout, camp::idx_t idx) {
    camp::idx_t v_idx = view_index(view, idx_seq_for(views));
    if(v_idx == -1) {
      std::cout << "View not found in tuple of views for this decision object\n";
      return false;
    }

    user_decided[idx][v_idx] = true;
    return true;
  }
  
};



template <camp::idx_t NumKernels, typename...Views>
struct LayoutDecisions {

  constexpr static camp::idx_t NumViews = sizeof...(Views);
  using LayoutTypes = camp::tuple<typename std::remove_reference<Views>::type::layout_type...>;
   
  std::array<LayoutTypes, NumKernels + 1> moments;
  camp::tuple<Views...> views;

  template <camp::idx_t...Is>
  LayoutDecisions(camp::tuple<Views...> _views, camp::idx_seq<Is...>) :views(_views) {
    
    std::cout << "instatiating a LayoutDecisions object with " << NumViews << " views\n";
  
    for(auto i = 0; i < moments.size(); i++) {
      camp::sink(add_decision(get<Is>(views), get<Is>(views).get_layout(), i)...);
    }

     
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
        //std::cout << "The views are refreences to the same views\n";
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
  bool add_decision(ViewType & view_to_change, LayoutType new_layout, camp::idx_t when) {
    update_moment(moments[when], view_to_change, new_layout, idx_seq_for(moments[when]));
    return 1;
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

  return DecisionType(views, idx_seq_for(views));
}


} //namespace RAJA 


#endif
