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

    isl_union_set * space = isl_union_set_read_from_str(ctx, s.str().c_str());

  return space;
}

template <camp::idx_t NumDims> 
isl_union_map * objective_function(isl_ctx * ctx, TermList<NumDims> * vars, Constraint<NumDims> * constraints, Expression<NumDims> * cost_expression) {

  auto obj = new NamedVariable<NumDims>("obj");

  auto obj_func = new EqualityConstraint<NumDims>(obj, cost_expression);

  auto out_vars = new TermList<NumDims>();
  out_vars->add(obj);
  for(auto var : vars->variables) {
    out_vars->add(var);
  }

  std::stringstream s; 
  //s << "{" << vars->rep() << " -> " << out_vars->rep() << " : " << constraints->rep() << " && " << obj_func->rep() << "}";
  s << "{" << vars->rep() << " -> " << out_vars->rep() << " : " << obj_func->rep() << "}";
  isl_union_map * map = isl_union_map_read_from_str(ctx, s.str().c_str());

  return map;
}






} //namespace RAJA 

#endif

