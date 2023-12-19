#ifndef _RAJA_LAYOUTISLTRANSLATION_HPP
#define _RAJA_LAYOUTISLTRANSLATION_HPP

#include "RAJA/loopchain/FormatDecisions/Model/ConstraintGrammar.hpp"
#include "RAJA/util/all-isl.hpp"
#include <sstream>

namespace RAJA {



template <camp::idx_t NumDims> 
isl_union_set * decision_space(isl_ctx * ctx, TermList<NumDims> * vars, Constraint<NumDims> * constraints) {


  
  std::stringstream s;
  s << "{ " << vars->rep() << " : " << constraints->rep() << " } ";
  std::cerr << "Decision space isl string: " << s.str() << "\n";
  isl_union_set * space = isl_union_set_read_from_str(ctx, s.str().c_str());


  return space;
}

template <camp::idx_t NumDims> 
isl_union_pw_qpolynomial * objective_function(isl_ctx * ctx, TermList<NumDims> * vars, Expression<NumDims> * cost_expression) {

  std::stringstream s; 
  //std::cout << "objective_function.\n" << std::flush;
  //std::cout << "variables: " << vars->rep() << "\n" << std::flush ;
  //std::cout << "costs: " << cost_expression->rep() << "\n" << std::flush ;
  s << "{" << vars->rep() << " -> " << cost_expression->rep() << "}";
  isl_union_pw_qpolynomial * map = isl_union_pw_qpolynomial_read_from_str(ctx, s.str().c_str());

  return map;
}





} //namespace RAJA 

#endif

