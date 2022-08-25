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
  std::cout << "\n\n" << s.str() << "\n";
  isl_union_set * space = isl_union_set_read_from_str(ctx, s.str().c_str());


  return space;
}

template <camp::idx_t NumDims> 
isl_union_pw_qpolynomial * objective_function(isl_ctx * ctx, TermList<NumDims> * vars, Expression<NumDims> * cost_expression) {

  std::stringstream s; 
  //s << "{" << vars->rep() << " -> " << out_vars->rep() << " : " << constraints->rep() << " && " << obj_func->rep() << "}";
  s << "{" << vars->rep() << " -> " << cost_expression->rep() << "}";
  isl_union_pw_qpolynomial * map = isl_union_pw_qpolynomial_read_from_str(ctx, s.str().c_str());

  return map;
}

isl_point * find_best_choice(isl_union_set * choices, isl_union_pw_qpolynomial * objective) {
  struct OptStruct {
    isl_union_pw_qpolynomial * obj;
    isl_point * lowest_point;
    camp::idx_t lowest_val;
  } fn_struct;
  fn_struct.obj = objective;
  fn_struct.lowest_point = NULL;
  fn_struct.lowest_val = INT_MAX;

  auto fn = [](isl_point * pnt, void * user) {
    OptStruct * fn_struct = (OptStruct*) user;
    isl_union_pw_qpolynomial * obj = isl_union_pw_qpolynomial_copy(fn_struct->obj);
    isl_point * lowest_point = fn_struct->lowest_point;
    camp::idx_t lowest_val = fn_struct->lowest_val;

    std::cout << "Applying polynomial.\n";
    isl_val * val = isl_union_pw_qpolynomial_eval(obj, isl_point_copy(pnt));
    std::cout << "Applied polynomial.\n";
    camp::idx_t v =  isl_val_get_num_si(val);
    std::cout << "val:" << v << "\n";
    if ( v < lowest_val) {
      fn_struct->lowest_val = v;
      fn_struct->lowest_point = pnt;
    }
    return isl_stat_ok;
  };

  isl_union_set_foreach_point(choices, fn, &fn_struct);

  return fn_struct.lowest_point;
}





} //namespace RAJA 

#endif

