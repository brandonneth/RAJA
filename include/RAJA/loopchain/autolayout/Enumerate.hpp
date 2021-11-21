#ifndef _RAJA_LCENUMERATE_HPP
#define _RAJA_LCENUMERATE_HPP

#include "RAJA/loopchain/SymExec.hpp"
#include <algorithm>
namespace RAJA {


template <typename ViewType>
auto enumerate_perms(ViewType v) {

  constexpr camp::idx_t NumDims = v.get_layout().n_dims;

  using ArrayType = std::array<camp::idx_t, NumDims>;

  std::vector<ArrayType> perms;

  ArrayType currPerm;

   
  for(int i = 0; i < NumDims; i++) {
    currPerm[i] = layout_to_perm(v.get_layout())[i];
  }    

    std::sort(currPerm.begin(), currPerm.end());

  do {
    ArrayType perm;
    for(int i = 0; i < NumDims; i++) {
      perm[i] = currPerm[i];
    }
    perms.push_back(perm);
  } while (std::next_permutation(currPerm.begin(), currPerm.end()));
  
  
  return perms;  
  
}
template <camp::idx_t CurrentMoment, camp::idx_t ViewIdx, typename DecisionType, typename...ViewTypes>
auto enumerate_decisions_for_views(std::vector<DecisionType> partial_decisions, tuple<ViewTypes...> views);

template <camp::idx_t CurrentMoment, camp::idx_t ViewIdx, typename DecisionType,  typename...ViewTypes>
auto enumerate_decisions_for_views(std::vector<DecisionType> partial_decisions, tuple<ViewTypes...> views) {
   if constexpr(ViewIdx == sizeof...(ViewTypes)) {
    
    return partial_decisions; 
  } else {  
    auto view_to_decide = get<ViewIdx>(views);

    decltype(view_to_decide) & view_ref = get<ViewIdx>(views);
    
    auto perm_choices = enumerate_perms(view_to_decide);
    std::array<camp::idx_t, view_to_decide.get_layout().n_dims> sizes;
    for(int i = 0; i < view_to_decide.get_layout().n_dims; i++) {
      sizes[i] = view_to_decide.get_layout().sizes[i];
    }



    std::vector<DecisionType> extended_decisions;
    for(auto &decision : partial_decisions) {
      for(auto perm : perm_choices) {
        DecisionType extended(decision);


        auto layout = make_permuted_layout(sizes, perm);
        extended.add_decision(view_ref, layout, CurrentMoment);
        extended_decisions.push_back(extended);
         
      }
    }
    //std::cout << "First decision, partially decided: " << extended_decisions[0] << "\n";
    return enumerate_decisions_for_views<CurrentMoment,ViewIdx+1>(extended_decisions, views);
  }

}

template <camp::idx_t MomentsDecided, typename DecisionType, typename...ViewTypes>
auto enumerate_moment(std::vector<DecisionType> partial_decisions, tuple<ViewTypes...> views) {
  std::vector<DecisionType> extended_decisions;
  
  for(DecisionType decision : partial_decisions) {
   
    std::vector<DecisionType> paths;
    paths.push_back(decision); 
    auto extended_paths = enumerate_decisions_for_views<MomentsDecided,0>(paths, views);
    for(auto & newPath : extended_paths) {
      extended_decisions.push_back(newPath);
    }
    //extended_decisions.insert(extended_decisions.end(), extended_paths.begin(), extended_paths.end());
  }
  return extended_decisions;
}


template <camp::idx_t MomentsDecided, camp::idx_t NumKernels, typename...ViewTypes>
auto enumerate_possible_decisions_helper(std::vector<LayoutDecisions<NumKernels,ViewTypes...>> partial_decisions, tuple<ViewTypes...> views) {
  //std::cout << "enumerate_possible_decisions_helper view is reference: " << std::is_reference<decltype(get<0>(views))>::value << "\n";

  if constexpr (MomentsDecided == NumKernels + 1 ) {
    return partial_decisions;
  } else {
    //add decisions for this moment
    auto more_decided = enumerate_moment<MomentsDecided>(partial_decisions, views );
    return enumerate_possible_decisions_helper<MomentsDecided+1, NumKernels>(more_decided, views);
  }
  
}

template <camp::idx_t NumKernels, typename... ViewTypes>
auto enumerate_possible_decisions(tuple<ViewTypes...> views) {

  std::cout << "enumerate_possible_decisions view is reference: " << std::is_reference<decltype(get<0>(views))>::value << "\n";
  
  using DecisionType = LayoutDecisions<NumKernels, ViewTypes...>;

  DecisionType no_decisions_made = create_layout_decisions<NumKernels>(views);

  std::vector<DecisionType> all_options;
  all_options.push_back(no_decisions_made);
  
  return enumerate_possible_decisions_helper<0, NumKernels>(all_options, views);
}






} //namespace RAJA

#endif


