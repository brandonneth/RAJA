#ifndef _RAJA_LAYOUTISLTRANSLATION_HPP
#define _RAJA_LAYOUTISLTRANSLATION_HPP

#include "RAJA/loopchain/autolayout/model/ConstraintGrammar.hpp"
#include "RAJA/loopchain/autolayout/model/DecisionVariable.hpp"
#include "RAJA/loopchain/autolayout/model/Constraint.hpp"
#include "RAJA/util/all-isl.hpp"


#include <sstream>

namespace RAJA {



template <camp::idx_t NumDims> 
isl_union_set * decision_space(isl_ctx * ctx, DecisionVariableList<NumDims> * vars, Constraint<NumDims> * constraints) {

  std::stringstream s;
  s << "{ " << vars->rep() << " : " << constraints->rep() << " } ";

  std::cout << "String for decision space\n";
  std::cout << s.str() << "\n";
  isl_union_set * space = isl_union_set_read_from_str(ctx, s.str().c_str());

  return space;
}







} //namespace RAJA 

#endif


