#ifndef _RAJA_MODELCONSTRAINT_HPP
#define _RAJA_MODELCONSTRAINT_HPP

#include <sstream>
#include "RAJA/loopchain/FormatDecisions/Enumerate.hpp"
#include "RAJA/loopchain/FormatDecisions/Model/ConstraintGrammar.hpp"
#include "RAJA/loopchain/FormatDecisions/Model/DecisionVariable.hpp"
namespace RAJA {
 

//Encodes constraint that for a particular computation, each view can only have one layout. 
//This is achieved by constraining the sum of the view's format variables for the kernel to 1. 
//This way only one format variable can be true and the others must be false.
template <camp::idx_t NumDims, typename Container>
EqualityConstraint<NumDims> * one_format_constraint(std::string viewName, camp::idx_t format_idx, Container possiblePerms) {
  std::vector<Expression<NumDims>*> format_vars;
  for(auto perm : possiblePerms) {
    std::vector<camp::idx_t> p;
    for(int i = 0; i < NumDims; i++) {
      p.push_back(perm[i]);
    }
    format_vars.push_back(new FormatDecision<NumDims>(viewName, format_idx, p));
  }

  auto sum = new SumExpression(format_vars);

  auto one = new Constant<NumDims>(1);
  return new EqualityConstraint<NumDims>(sum, one);  
  

} 

//Encodes constraint that at each conversion point, only one conversion can happen.
//Similar to the one_format_constraint, we constraint the sum to 1.
template <camp::idx_t NumDims>
EqualityConstraint<NumDims> * one_conversion_constraint(std::string viewName, camp::idx_t conv_idx) {
      
  
  if(conv_idx == 0) {
    std::cerr << "one_conversion_constraint conv_idx is 0, which is not allowed. "
                 "The first index on conversions is 1.\n";
    return new EqualityConstraint<NumDims>(new Constant<NumDims>(0), new Constant<NumDims>(0));
  }    
  std::vector<Expression<NumDims>*> conversion_vars;
  auto all_perms = enumerate_perms<NumDims>();
  for(auto perm1 : all_perms) {
    for(auto perm2 : all_perms) {
      std::vector<camp::idx_t> p1;
      std::vector<camp::idx_t> p2;
      for(int i = 0; i < NumDims; i++) {
        p1.push_back(perm1[i]);
        p2.push_back(perm2[i]);
      }   
      conversion_vars.push_back(new ConversionDecision<NumDims>(viewName, conv_idx, p1, p2));
    }
  }

  auto sum = new SumExpression(conversion_vars);

  auto one = new Constant<NumDims>(1);
  return new EqualityConstraint<NumDims>(sum, one);  
  

}

//Encodes constraint that for the conversion from one format to the next that the
// original format matches the input format of the conversion stage. This is done
// by setting the sum of the conversion decisions with a particular input to the
// value of the format decision for that input format (conversion2.in == format1). 
template <camp::idx_t NumDims>
AndConstraint<NumDims>* conversion_input_constraint(std::string viewName, camp::idx_t conv_idx) {
 
  AndConstraint<NumDims>* constraints = new AndConstraint<NumDims>();
  if(conv_idx == 0) {
    std::cerr << "conversion_input_constraint passed a conversion index of 0. The first conversion index is 1.\n";
    auto eq = new EqualityConstraint<NumDims>(new Constant<NumDims>(3), new Constant<NumDims>(3));
    constraints->add(eq); 
    return constraints;
  }
  auto all_perms = enumerate_perms<NumDims>();
      
  for(auto perm : all_perms) {
    std::vector<camp::idx_t> in_p;
    for(int i = 0; i < NumDims; i++) {
      in_p.push_back(perm[i]);
    }

    FormatDecision<NumDims> * format = new FormatDecision<NumDims>(viewName, conv_idx-1, in_p);

    auto sum = new SumExpression<NumDims>();

    for(auto out_perm : all_perms) {
      std::vector<camp::idx_t> out_p;
      for(int i = 0; i < NumDims; i++) {out_p.push_back(out_perm[i]);}
      
      sum->add(new ConversionDecision<2>(viewName, conv_idx, in_p, out_p));
    }
    
    constraints->add(new EqualityConstraint<NumDims>(format, sum));           
  }

  return constraints;
      
}

//Encodes complementary constraint to conversion_input_constraint, that the output format of a conversion matches the format decision for the moment after the conversion. 
template <camp::idx_t NumDims>
AndConstraint<NumDims>* conversion_output_constraint(std::string viewName, camp::idx_t conv_idx) {

  if(conv_idx == 0) {
    std::cerr << "conversion_input_constraint passed a conversion index of 0. The first conversion index is 1.\n";
    auto c =  new EqualityConstraint<NumDims>(new Constant<NumDims>(0), new Constant<NumDims>(0));
    return new AndConstraint<NumDims>(c,c);
  }    
  auto all_perms = enumerate_perms<NumDims>();
      
  AndConstraint<NumDims>* constraints = new AndConstraint<NumDims>();
  for(auto perm : all_perms) {
    std::vector<camp::idx_t> out_p;
    for(int i = 0; i < NumDims; i++) {
      out_p.push_back(perm[i]);
    }

    FormatDecision<NumDims> * format = new FormatDecision<NumDims> (viewName, conv_idx, out_p);

    auto sum = new SumExpression<NumDims>();

    for(auto in_perm : all_perms) {
      std::vector<camp::idx_t> in_p;
      for(int i = 0; i < NumDims; i++) {in_p.push_back(in_perm[i]);}
      
      sum->add(new ConversionDecision<2>(viewName, conv_idx, in_p, out_p));
    }
    
    constraints->add(new EqualityConstraint<NumDims>(format, sum));           
  }

  return constraints;
      
}

//Encodes the truth-value natuer of the decision variables. Each variable has to be 1 or 0.
template <camp::idx_t NumDims, typename Container>
AndConstraint<NumDims>* binary_choice_constraints(std::string viewName, camp::idx_t idx, Container possiblePerms) {
  AndConstraint<NumDims>* constraints = new AndConstraint<NumDims>();
  
  std::vector<Expression<NumDims>*> vars;
  for(auto perm : possiblePerms) {
    std::vector<camp::idx_t> p;
    for(int i = 0; i < NumDims; i++) {
      p.push_back(perm[i]);
    }
    vars.push_back(new FormatDecision<NumDims>(viewName, idx, p));
  }

  auto zero = new Constant<NumDims>(0);
  auto one = new Constant<NumDims>(1);
  for(auto v : vars) {
   constraints->add(new LeqConstraint<NumDims>(v, one));
   constraints->add(new GeqConstraint<NumDims>(v, zero));
  }
  return constraints;   
}


//Collects encodes all the constraints above for all the decision variables for a particular number of kernels.
template <camp::idx_t NumDims, typename Container>
AndConstraint<NumDims>* base_constraints(std::string viewName, camp::idx_t num_knls, Container possiblePerms) {
  AndConstraint<NumDims>* constraints = new AndConstraint<NumDims>();
  
  constraints->add(one_format_constraint<NumDims>(viewName, 0, possiblePerms));
  constraints->add(binary_choice_constraints<NumDims>(viewName, 0, possiblePerms));

  for(auto i = 1; i <= num_knls+1; i++) {
    constraints->add(one_format_constraint<NumDims>(viewName, i, possiblePerms));
    //constraints->add(one_conversion_constraint<NumDims>(viewName, i));
    //constraints->add(conversion_input_constraint<NumDims>(viewName, i));
    //constraints->add(conversion_output_constraint<NumDims>(viewName, i));
    constraints->add(binary_choice_constraints<NumDims>(viewName, i, possiblePerms));
  }
  return constraints;
}
} //namespace RAJA

#endif

