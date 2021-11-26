#ifndef _RAJA_MODELCONSTRAINTGRAMMAR_HPP
#define _RAJA_MODELCONSTRAINTGRAMMAR_HPP

#include <sstream>

#include "RAJA/loopchain/autolayout/Enumerate.hpp"

namespace RAJA {

  template <camp::idx_t I>
  class Term {
    public:
    virtual std::string rep() const;
  };

  template <camp::idx_t I>
  class Constraint {
    public:
    virtual std::string rep();
  };

  template <camp::idx_t I>
  class EqualityConstraint : public Constraint<I> {
    public:
    Term<I> * lhs;
    Term<I> * rhs;

    EqualityConstraint(Term<I> * lhs_, Term<I> * rhs_) : lhs(lhs_), rhs(rhs_) {}

    virtual std::string rep() {
      std::stringstream s;
      s << lhs->rep() << " = " << rhs->rep();
      return s.str();
    }
  };

  template <camp::idx_t I>
  class LeqConstraint : public Constraint<I> {
    public:
    Term<I> * lhs;
    Term<I> * rhs;

    LeqConstraint(Term<I> * lhs_, Term<I> * rhs_) : lhs(lhs_), rhs(rhs_) {}

    virtual std::string rep() {
      std::stringstream s;
      s << lhs->rep() << " <= " << rhs->rep();
      return s.str();
    }
  };

  template <camp::idx_t I>
  class GeqConstraint : public Constraint<I> {
    public:
    Term<I> * lhs;
    Term<I> * rhs;

    GeqConstraint(Term<I> * lhs_, Term<I> * rhs_) : lhs(lhs_), rhs(rhs_) {}

    virtual std::string rep() {
      std::stringstream s;
      s << lhs->rep() << " >= " << rhs->rep();
      return s.str();
    }
  };

  template <camp::idx_t I>
  class AndConstraint : public Constraint<I> {
    public:
    std::vector<Constraint<I>*> constraints;

  
    AndConstraint() = default;
    AndConstraint(Constraint<I>* c1, Constraint<I>* c2) {
      constraints.push_back(c1);
      constraints.push_back(c2);
    }
    AndConstraint(std::vector<Constraint<I>*> cs) : constraints(cs) {}


    void add(Constraint<I>* c) {
      constraints.push_back(c);
    }
    virtual std::string rep() {
      if(constraints.size() == 0) {
        return "";
      }
      if (constraints.size() == 1) {
        return constraints[0]->rep();
      }

      std::stringstream s;
      s << constraints[0]->rep();
      for(int i = 1; i < constraints.size(); i++) {
        s << " && " << constraints[i]->rep();
      }
      return s.str();

    }
  };


  template <camp::idx_t I>
  class DecisionVariable : public Term<I> {
    public:
        
    virtual std::string rep() const {
       return "DECISION VARIABLE REP FUNCTION BEING CALLED IN ERROR. You're probably using an object instead of a pointer.";;
    }
   };
 
  
   
  template <camp::idx_t I>
  class DecisionVariableSum : public Term<I> {
    public:
    std::vector<DecisionVariable<I>*> variables;

    DecisionVariableSum() = default;
    DecisionVariableSum(std::vector<DecisionVariable<I>*> v_) : variables(v_) {}
    void add(DecisionVariable<I>* v) {
      variables.push_back(v);
    }

    virtual std::string rep() const {
      std::stringstream s;
      s << variables[0]->rep();
      for(int i = 1; i < variables.size(); i++) {
        s << " + " << variables[i]->rep();
      }
      return s.str();
    }

  };

  template <camp::idx_t I>
  class DecisionVariableList : public Term<I> {
    public:
    std::vector<DecisionVariable<I>*> variables;

    DecisionVariableList() = default;
    DecisionVariableList(std::vector<DecisionVariable<I>*> v_) : variables(v_) {}
    
    void add(DecisionVariable<I>* v) {
      variables.push_back(v);
    }

    virtual std::string rep() const {
      std::stringstream s;
      s << "[ ";
      s << variables[0]->rep();
      for(int i = 1; i < variables.size(); i++) {
        s << " , " << variables[i]->rep();
      }
      s << " ]";
      return s.str();
    }

  };

  template <camp::idx_t I>
  class Constant : public Term<I> {
    public:
    camp::idx_t c;
    Constant(camp::idx_t c_) : c(c_) {}

    virtual std::string rep() const {
      std::stringstream s;
      s << c;
      return s.str();
    }
  };


  


} //namespace RAJA

#endif


