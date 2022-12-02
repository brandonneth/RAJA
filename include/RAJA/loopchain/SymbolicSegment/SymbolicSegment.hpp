#ifndef RAJA_SymbolicSegment_HPP
#define RAJA_SymbolicSegment_HPP

#include "RAJA/config.hpp"
#include "RAJA/util/Operators.hpp"
#include "RAJA/util/Layout.hpp"
#include "RAJA/loopchain/PolicyToNestingOrder.hpp"

namespace RAJA 
{

struct SymbolicSegment;
struct SymbolicBound;

struct SymbolicBound {
  
  SymbolicBound() {}

  virtual idx_t const value() {
    return 0;
  }
  
}; //SymbolicBound

struct SymbolicBoundConstant :SymbolicBound {
  idx_t c;

  SymbolicBoundConstant(idx_t val) : c(val) {}

  virtual idx_t const value() {
    return c;
  }
};

struct SymbolicBoundSegment : SymbolicBound {
  const SymbolicSegment & seg;

  SymbolicBoundSegment(const SymbolicSegment &s) : seg(s) {}

  virtual idx_t const value();  
};
struct SymbolicSegment {
 
  SymbolicBound *low;
  SymbolicBound *high; 
 
  idx_t currVal;
  SymbolicSegment() : low(NULL), high(NULL) {
  }

  SymbolicSegment(SymbolicBound * l, SymbolicBound * h) : low(l), high(h) {
    currVal = low->value();
  }
  SymbolicBound * operator * (idx_t c);

  idx_t value() const {
    return currVal;
  }
}; //SymbolicSegment




struct SymbolicBoundAdd {
  SymbolicBound *left;
  SymbolicBound * right;

  idx_t value() {
    return left->value() + right->value();
  }
};
struct SymbolicBoundMult : SymbolicBound {
  SymbolicBound * left;
  SymbolicBound * right;

  SymbolicBoundMult(SymbolicBound * l, SymbolicBound * r) : left(l), right(r) {}
  virtual idx_t const value() {
    return left->value() * right->value();
  }
};


idx_t const SymbolicBoundSegment::value() {
    return seg.value();
}

SymbolicBound * SymbolicSegment::operator * (idx_t c) {
  return new SymbolicBoundMult(new SymbolicBoundSegment(*this), new SymbolicBoundConstant(c));
}



SymbolicBound * make_symbolic_bound(idx_t x) {
  return new SymbolicBoundConstant(x);
}

SymbolicBound * make_symbolic_bound(int x) {
  return new SymbolicBoundConstant(x);
}

SymbolicBound * make_symbolic_bound(const SymbolicSegment &x) {
  return new SymbolicBoundSegment(x);
}
SymbolicBound * make_symbolic_bound(SymbolicBound * x) {
  return x;
}
template <typename T1, typename T2>
SymbolicSegment make_symbolic_segment(T1 l, T2 r) {
  
  auto low = make_symbolic_bound(l);
  auto high = make_symbolic_bound(r); 
  return SymbolicSegment(low,high);
}


} //namespace RAJA

#endif 
