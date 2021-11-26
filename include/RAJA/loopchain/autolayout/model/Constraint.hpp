#ifndef _RAJA_MODELCONSTRAINT_HPP
#define _RAJA_MODELCONSTRAINT_HPP

#include <sstream>

#include "RAJA/loopchain/autolayout/model/ConstraintGrammar.hpp"
#include "RAJA/loopchain/autolayout/model/DecisionVariable.hpp"
namespace RAJA {
 

template <camp::idx_t NumDims>
EqualityConstraint<NumDims> * one_format_constraint(std::string viewName, camp::idx_t idx) {
  std::vector<DecisionVariable<NumDims>*> format_vars;
  auto all_perms = enumerate_perms<NumDims>();
  for(auto perm : all_perms) {
    std::vector<camp::idx_t> p;
    for(int i = 0; i < NumDims; i++) {
      p.push_back(perm[i]);
    }
    format_vars.push_back(new FormatDecision<NumDims>(viewName, idx, p));
  }

  auto sum = new DecisionVariableSum(format_vars);

  auto one = new Constant<NumDims>(1);
  return new EqualityConstraint<NumDims>(sum, one);  
  

} 

template <camp::idx_t NumDims>
EqualityConstraint<NumDims> * one_conversion_constraint(std::string viewName, camp::idx_t idx) {
      
  if(idx == 0) {
    return new EqualityConstraint<NumDims>(new Constant<NumDims>(0), new Constant<NumDims>(0));
  }    
  std::vector<DecisionVariable<NumDims>*> conversion_vars;
  auto all_perms = enumerate_perms<NumDims>();
  for(auto perm1 : all_perms) {
    for(auto perm2 : all_perms) {
      std::vector<camp::idx_t> p1;
      std::vector<camp::idx_t> p2;
      for(int i = 0; i < NumDims; i++) {
        p1.push_back(perm1[i]);
        p2.push_back(perm2[i]);
      }   
      conversion_vars.push_back(new ConversionDecision<NumDims>(viewName, idx, p1, p2));
    }
  }

  auto sum = new DecisionVariableSum(conversion_vars);

  auto one = new Constant<NumDims>(1);
  return new EqualityConstraint<NumDims>(sum, one);  
  

}

template <camp::idx_t NumDims>
AndConstraint<NumDims>* conversion_input_constraint(std::string viewName, camp::idx_t idx) {
 
  auto all_perms = enumerate_perms<NumDims>();
      
  AndConstraint<NumDims>* constraints = new AndConstraint<NumDims>();
  for(auto perm : all_perms) {
    std::vector<camp::idx_t> p;
    for(int i = 0; i < NumDims; i++) {
      p.push_back(perm[i]);
    }

    FormatDecision<NumDims> * format = new FormatDecision<NumDims> (viewName, idx, p);

    auto sum = new DecisionVariableSum<NumDims>();

    for(auto out_perm : all_perms) {
      std::vector<camp::idx_t> out_p;
      for(int i = 0; i < NumDims; i++) {out_p.push_back(out_perm[i]);}
      
      sum->add(new ConversionDecision<2>(viewName, idx+1, p, out_p));
    }
    
    constraints->add(new EqualityConstraint<NumDims>(format, sum));           
  }

  return constraints;
      
}
template <camp::idx_t NumDims>
AndConstraint<NumDims>* conversion_output_constraint(std::string viewName, camp::idx_t idx) {

  if(idx == 0) {
    auto c =  new EqualityConstraint<NumDims>(new Constant<NumDims>(0), new Constant<NumDims>(0));
    return new AndConstraint<NumDims>(c,c);
  }    
  auto all_perms = enumerate_perms<NumDims>();
      
  AndConstraint<NumDims>* constraints = new AndConstraint<NumDims>();
  for(auto perm : all_perms) {
    std::vector<camp::idx_t> p;
    for(int i = 0; i < NumDims; i++) {
      p.push_back(perm[i]);
    }

    FormatDecision<NumDims> * format = new FormatDecision<NumDims> (viewName, idx, p);

    auto sum = new DecisionVariableSum<NumDims>();

    for(auto in_perm : all_perms) {
      std::vector<camp::idx_t> in_p;
      for(int i = 0; i < NumDims; i++) {in_p.push_back(in_perm[i]);}
      
      sum->add(new ConversionDecision<2>(viewName, idx, in_p, p));
    }
    
    constraints->add(new EqualityConstraint<NumDims>(format, sum));           
  }

  return constraints;
      
}

template <camp::idx_t NumDims>
AndConstraint<NumDims>* binary_choice_constraints(std::string viewName, camp::idx_t idx) {
  AndConstraint<NumDims>* constraints = new AndConstraint<NumDims>();
  
  std::vector<DecisionVariable<NumDims>*> vars;
  auto all_perms = enumerate_perms<NumDims>();
  for(auto perm : all_perms) {
    std::vector<camp::idx_t> p;
    for(int i = 0; i < NumDims; i++) {
      p.push_back(perm[i]);
    }
    vars.push_back(new FormatDecision<NumDims>(viewName, idx, p));
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
        vars.push_back(new ConversionDecision<NumDims>(viewName, idx, p1, p2));
      }
    }
  }

  auto zero = new Constant<2>(0);
  auto one = new Constant<2>(1);
  for(auto v : vars) {
   constraints->add(new LeqConstraint<NumDims>(v, one));
   constraints->add(new GeqConstraint<NumDims>(v, zero));
  }



  return constraints;
      
}


template <camp::idx_t NumDims>
AndConstraint<NumDims>* all_constraints(std::string viewName, camp::idx_t num_knls) {
  AndConstraint<NumDims>* constraints = new AndConstraint<NumDims>();
  
    constraints->add(one_format_constraint<NumDims>(viewName, 0));
    constraints->add(binary_choice_constraints<NumDims>(viewName, 0));
    constraints->add(conversion_input_constraint<NumDims>(viewName, 0));

  for(auto i = 1; i < num_knls; i++) {
    constraints->add(one_format_constraint<NumDims>(viewName, i));
    constraints->add(one_conversion_constraint<NumDims>(viewName, i));
    constraints->add(conversion_input_constraint<NumDims>(viewName, i));
    constraints->add(conversion_output_constraint<NumDims>(viewName, i));
    constraints->add(binary_choice_constraints<NumDims>(viewName, i));
  }
    constraints->add(one_format_constraint<NumDims>(viewName, num_knls));
    constraints->add(one_conversion_constraint<NumDims>(viewName, num_knls));
  
    constraints->add(conversion_output_constraint<NumDims>(viewName, num_knls));
    constraints->add(binary_choice_constraints<NumDims>(viewName, num_knls));
  return constraints;
}
} //namespace RAJA

#endif


