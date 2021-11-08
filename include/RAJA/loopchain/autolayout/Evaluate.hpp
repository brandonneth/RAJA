#ifndef _RAJA_LCEVALUATE_HPP
#define _RAJA_LCEVALUATE_HPP

#include "RAJA/loopchain/SymExec.hpp"
#include <ctime>
namespace RAJA {


template <typename DecisionType, typename Knl1, typename Knl2>
auto evaluate_decisions(std::vector<DecisionType> decision_objects, Knl1 knl1, Knl2 knl2) {
  
  auto comp1 = manipulate_layouts(decision_objects[0], knl1, knl2);
  for(auto decision : decision_objects) {
    
    auto computation = manipulate_layouts(decision, knl1, knl2);

    auto start = std::clock();
    computation();
    auto stop = std::clock();

    auto t = stop - start;

    std::cout << decision << ", " << t << "\n";
    
  }  

  return comp1;
}





} //namespace RAJA

#endif


