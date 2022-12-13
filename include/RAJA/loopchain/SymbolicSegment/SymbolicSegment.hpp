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

  virtual idx_t value() const {
    return 0;
  }
  
}; //SymbolicBound

struct SymbolicBoundConstant :SymbolicBound {
  idx_t c;

  SymbolicBoundConstant(idx_t val) : c(val) {}

  virtual idx_t  value() const {
    return c;
  }
};

struct SymbolicBoundSegment : SymbolicBound {
  const SymbolicSegment & seg;

  SymbolicBoundSegment(const SymbolicSegment &s) : seg(s) {}

  virtual idx_t value() const;  
};


struct SymbolicIter {
  SymbolicBound * bound;
  const SymbolicSegment * seg;
  Index_type offset;

  SymbolicIter(SymbolicBound * b, const SymbolicSegment * s) : bound(b), seg(s), offset(0) {}
  Index_type operator* () {
    return bound->value() + offset;
  }

  SymbolicIter operator + (Index_type i) {
    SymbolicIter newIter = SymbolicIter(bound, seg);
    newIter.offset = offset + i;
    return newIter;
  }
 
  Index_type operator - (SymbolicIter other) {
    return *(*this) - *(other);
  } 
  
  Index_type operator[](Index_type i) {
    return *(*this) + i;
  }

  void update(Index_type i);
};

struct SymbolicSegment {
  using iterator = SymbolicIter;
  using value_type = Index_type;
  using IndexType = Index_type;
  SymbolicBound *low;
  SymbolicBound *high; 
  std::string name;
    SymbolicSegment() : low(NULL), high(NULL) {
  }

  SymbolicSegment(SymbolicBound * l, SymbolicBound * h, std::string n) : low(l), high(h), name(n) {
    currVal = new IndexType[1];
    currVal[0] = low->value();
  }
  
  SymbolicBound* operator * (idx_t c);

  void update_offset(IndexType i) const {
    currVal[0] = low->value() + i;
  }
  idx_t value() const {
    return currVal[0];
  }

  SymbolicIter begin() const {
    return SymbolicIter(low, this);
  }

  SymbolicIter end() const {
    return SymbolicIter(high, this);
  }
  
  void increment() {
    currVal[0] += 1;
  }
  void decrement() {
    currVal[0] -= 1;
  }

  void reset() {
    currVal[0] = low->value();
  }
private:
idx_t * currVal;

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
  virtual idx_t value() const {
    return left->value() * right->value();
  }
};


idx_t SymbolicBoundSegment::value() const {
    return seg.value();
}

SymbolicBound * SymbolicSegment::operator * (idx_t c) {
  return new SymbolicBoundMult(new SymbolicBoundSegment(*this), new SymbolicBoundConstant(c));
}


void SymbolicIter::update(Index_type i) {
    seg->update_offset(i);
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
SymbolicSegment make_symbolic_segment(const T1 & l, const T2 r, std::string n) {
  
  SymbolicBound * low = make_symbolic_bound(l);
  SymbolicBound * high = make_symbolic_bound(r); 
  return SymbolicSegment(low,high,n);
}

template <typename T1, typename T2>
SymbolicSegment make_symbolic_segment(const T1 & l, const T2 r) {
  
  SymbolicBound * low = make_symbolic_bound(l);
  SymbolicBound * high = make_symbolic_bound(r); 
  return SymbolicSegment(low,high,"");
}
} //namespace RAJA

template <>
struct std::iterator_traits<RAJA::SymbolicIter> {
  typedef RAJA::Index_type value_type;
  typedef RAJA::Index_type difference_type;
  typedef std::random_access_iterator_tag iterator_category;
  typedef RAJA::Index_type* pointer;
  typedef RAJA::Index_type& reference;
};


#endif 
