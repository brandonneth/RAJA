#ifndef _RAJA_ModelConstraintGrammar_HPP
#define _RAJA_ModelConstraintGrammar_HPP

namespace RAJA
{


template <camp::idx_t I>
class Constraint
{
public:
  virtual std::string rep();
};


template <camp::idx_t I>
class Expression
{
public:
  virtual std::string rep() const;
};


template <camp::idx_t I>
class EqualityConstraint : public Constraint<I>
{
public:
  Expression<I>* lhs;
  Expression<I>* rhs;

  EqualityConstraint(Expression<I>* lhs_, Expression<I>* rhs_)
      : lhs(lhs_), rhs(rhs_)
  {
  }

  virtual std::string rep()
  {
    std::stringstream s;
    s << lhs->rep() << " = " << rhs->rep();
    return s.str();
  }
};

template <camp::idx_t I>
class LeqConstraint : public Constraint<I>
{
public:
  Expression<I>* lhs;
  Expression<I>* rhs;

  LeqConstraint(Expression<I>* lhs_, Expression<I>* rhs_) : lhs(lhs_), rhs(rhs_)
  {
  }

  virtual std::string rep()
  {
    std::stringstream s;
    s << lhs->rep() << " <= " << rhs->rep();
    return s.str();
  }
};

template <camp::idx_t I>
class GeqConstraint : public Constraint<I>
{
public:
  Expression<I>* lhs;
  Expression<I>* rhs;

  GeqConstraint(Expression<I>* lhs_, Expression<I>* rhs_) : lhs(lhs_), rhs(rhs_)
  {
  }

  virtual std::string rep()
  {
    std::stringstream s;
    s << lhs->rep() << " >= " << rhs->rep();
    return s.str();
  }
};

template <camp::idx_t I>
class AndConstraint : public Constraint<I>
{
public:
  std::vector<Constraint<I>*> constraints;


  AndConstraint() = default;
  AndConstraint(Constraint<I>* c1, Constraint<I>* c2)
  {
    constraints.push_back(c1);
    constraints.push_back(c2);
  }
  AndConstraint(std::vector<Constraint<I>*> cs) : constraints(cs) {}


  void add(Constraint<I>* c) { constraints.push_back(c); }
  virtual std::string rep()
  {
    if (constraints.size() == 0) {
      return "";
    }
    if (constraints.size() == 1) {
      return constraints[0]->rep();
    }

    std::stringstream s;
    s << constraints[0]->rep();
    for (size_t i = 1; i < constraints.size(); i++) {
      s << " && " << constraints[i]->rep();
    }
    return s.str();
  }
};


template <camp::idx_t I>
class Term : public Expression<I>
{
public:
  virtual std::string rep() const;
  virtual camp::idx_t estimate_cost();
};

template <camp::idx_t I>
class MultExpression : public Expression<I>
{
public:
  Expression<I>* lhs;
  Expression<I>* rhs;

  MultExpression(Expression<I>* l, Expression<I>* r) : lhs(l), rhs(r) {}
  virtual std::string rep() const
  {
    std::stringstream s;
    s << lhs->rep() << " * " << rhs->rep();
    return s.str();
  }
};

template <camp::idx_t I>
class Constant : public Expression<I>
{
public:
  camp::idx_t c;
  Constant(camp::idx_t c_) : c(c_) {}

  virtual std::string rep() const
  {
    std::stringstream s;
    s << c;
    return s.str();
  }
};

template <camp::idx_t I>
class NamedVariable : public Term<I>
{
public:
  std::string name;
  NamedVariable<I>(std::string n) : name(n) {}

  virtual std::string rep() const { return name; }
  virtual camp::idx_t estimate_cost() { return 0; }
};


template <camp::idx_t I>
class DecisionVariable : public Term<I>
{
public:
  virtual std::string rep() const
  {
    return "DECISION VARIABLE REP FUNCTION BEING CALLED IN ERROR. You're "
           "probably using an object instead of a pointer.";
    ;
  }


  virtual camp::idx_t estimate_cost() { return 0; }
};


template <camp::idx_t I>
class SumExpression : public Expression<I>
{
public:
  std::vector<Expression<I>*> variables;

  SumExpression() = default;
  SumExpression(std::vector<Expression<I>*> v_) : variables(v_) {}
  void add(Expression<I>* v) { variables.push_back(v); }

  virtual std::string rep() const
  {
    std::stringstream s;
    s << variables[0]->rep();
    for (size_t i = 1; i < variables.size(); i++) {
      s << " + " << variables[i]->rep();
    }
    return s.str();
  }
};

template <camp::idx_t I>
class TermList
{
public:
  std::vector<Term<I>*> variables;

  TermList() = default;
  TermList(std::vector<Term<I>*> v_) : variables(v_) {}

  void add(Term<I>* v) { variables.push_back(v); }

  void prepend(Term<I>* v) { variables.insert(variables.begin(), v); }

  virtual std::string rep() const
  {
    std::stringstream s;
    s << "[ ";
    s << variables[0]->rep();
    for (size_t i = 1; i < variables.size(); i++) {
      s << " , " << variables[i]->rep();
    }
    s << " ]";
    return s.str();
  }
};
}  // namespace RAJA

#endif