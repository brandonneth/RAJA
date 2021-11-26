#ifndef _RAJA_MODELDECISIONVARIABLE_HPP
#define _RAJA_MODELDECISIONVARIABLE_HPP

#include <sstream>
#include "RAJA/loopchain/autolayout/model/ConstraintGrammar.hpp"


namespace RAJA {

  
  template <camp::idx_t I>
  std::size_t estimate_cost(DecisionVariable<I> * d); 

  template <camp::idx_t I>
  class FormatDecision : public DecisionVariable<I> {
    public:
    std::string viewName;
    camp::idx_t idx;

    std::array<camp::idx_t, I> permutation;
    FormatDecision(std::string viewName_, camp::idx_t idx_, std::vector<camp::idx_t> perm) 
      : viewName(viewName_), idx(idx_)
    {
      for(auto i = 0; i < I; i++) {
        permutation[i] = perm.at(i);
      }
    }

    std::string rep() const {
      std::stringstream s;
      s << viewName << "_f" << idx;
      for(auto i = 0; i < permutation.size(); i++) {
        s << "_" << permutation[i];
      }
      return s.str();

    }

  };
    
  template <camp::idx_t I>
  class ConversionDecision : public DecisionVariable<I> {
   public:

    ConversionDecision(std::string viewName_, camp::idx_t idx_, std::vector<camp::idx_t> in_perm_, std::vector<camp::idx_t> out_perm_) : 
      viewName(viewName_), 
      idx(idx_) 
    {
      
 
    for(auto i = 0; i < I; i++) {
        in_perm[i] = in_perm_.at(i);
        out_perm[i] = out_perm_.at(i);
      }

    }

    std::string viewName;
    camp::idx_t idx;
    std::array<camp::idx_t,I> in_perm;
    std::array<camp::idx_t,I> out_perm;


    virtual std::string rep() const {

      std::stringstream s;
      s << viewName << "_conv" << idx << "_in";
      for(auto i = 0; i < in_perm.size(); i++) {
        s << "_" << in_perm[i];
      }
      s << "_out";
      for(auto i = 0; i < out_perm.size(); i++) {
        s << "_" << out_perm[i];
      }
      return s.str();
    }


  };


template <camp::idx_t NumDims>
DecisionVariableList<NumDims> * decision_variables(std::string viewName, camp::idx_t idx) {
  auto vars = new DecisionVariableList<NumDims>();
  auto all_perms = enumerate_perms<NumDims>();
  for(auto perm : all_perms) {
    std::vector<camp::idx_t> p;
    for(int i = 0; i < NumDims; i++) {
      p.push_back(perm[i]);
    }
    vars->add(new FormatDecision<NumDims>(viewName, idx, p));
  }

  //format only for idx=0
  if(idx != 0) {
    for(auto perm1 : all_perms) {
      for(auto perm2 : all_perms) {
        std::vector<camp::idx_t> p1;
        std::vector<camp::idx_t> p2;
        for(int i = 0; i < NumDims; i++) {
          p1.push_back(perm1[i]);
          p2.push_back(perm2[i]);
        }   
        vars->add(new ConversionDecision<NumDims>(viewName, idx, p1, p2));
      }
    }
  }

 


  return vars;
      


}

template <camp::idx_t NumDims>
DecisionVariableList<NumDims> * all_decision_variables(std::string viewName, camp::idx_t max_idx) {
  auto vars = new DecisionVariableList<NumDims>();
  auto all_perms = enumerate_perms<NumDims>();
 
for(int idx = 0; idx <= max_idx; idx++) {

  for(auto perm : all_perms) {
    std::vector<camp::idx_t> p;
    for(int i = 0; i < NumDims; i++) {
      p.push_back(perm[i]);
    }
    vars->add(new FormatDecision<NumDims>(viewName, idx, p));
  }

  //format only for idx=0
  if(idx != 0) {
    for(auto perm1 : all_perms) {
      for(auto perm2 : all_perms) {
        std::vector<camp::idx_t> p1;
        std::vector<camp::idx_t> p2;
        for(int i = 0; i < NumDims; i++) {
          p1.push_back(perm1[i]);
          p2.push_back(perm2[i]);
        }   
        vars->add(new ConversionDecision<NumDims>(viewName, idx, p1, p2));
      }
    }
  }
}
  return vars;  
}

#include <ctime>
template <camp::idx_t I>
std::size_t estimate_cost(ConversionDecision<I> d) {
  using V =  View<double, Layout<I>>;

  std::size_t n = 1024;
  //V in(d.in_perm);
  //V out(d.out_perm);

  /*
 * create and then time the kernel that does a small version of the conversion
 */

  return 0;

} 
} //namespace RAJA

#endif


