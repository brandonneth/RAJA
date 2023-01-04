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

  virtual std::ostream & write(std::ostream& s) {
    return s;
  }
  friend std::ostream& operator<<(std::ostream & s, SymbolicBound& b) {
    return b.write(s);
  }
  
}; //SymbolicBound

struct SymbolicBoundConstant :SymbolicBound {
  idx_t c;

  SymbolicBoundConstant(idx_t val) : c(val) {}

  idx_t  value() const {
    return c;
  }

  std::ostream& write(std::ostream & s) {
    s << c;
    return s;
  }
};

struct SymbolicBoundSegment : SymbolicBound {
  const SymbolicSegment & seg;

  SymbolicBoundSegment(const SymbolicSegment &s) : seg(s) {}

  virtual idx_t value() const;  

  std::ostream& write(std::ostream & s); 
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
  size_t numThreads;
  SymbolicSegment() : low(NULL), high(NULL) {}

  SymbolicSegment(SymbolicBound * l, SymbolicBound * h) : low(l), high(h), name(std::string()) {
    init_val();
  }

  SymbolicSegment(SymbolicBound * l, SymbolicBound * h, std::string n) : low(l), high(h), name(n) {
    init_val();
      }

  void init_val() {
#ifdef RAJA_ENABLE_OPENMP
    numThreads = omp_get_max_threads();
#else
    numThreads = 1;
#endif
    currVal = new IndexType[numThreads];
    for(size_t i = 0; i < numThreads; i++) {
      currVal[i] = low->value();
    }
 
    
  }
  
  SymbolicBound* operator * (idx_t c);
  SymbolicBound* operator + (idx_t c);

  void update_offset(IndexType i) const {
#ifdef RAJA_ENABLE_OPENMP
    int thread = omp_get_thread_num();
#else
    int thread = 0;
#endif
   currVal[thread] = low->value() + i;
  }

  idx_t value() const {
#ifdef RAJA_ENABLE_OPENMP
    int thread = omp_get_thread_num();
#else
    int thread = 0;
#endif
   return currVal[thread];
  }

  SymbolicIter begin() const {
    return SymbolicIter(low, this);
  }

  SymbolicIter end() const {
    return SymbolicIter(high, this);
  }
  
  void increment() {
#ifdef RAJA_ENABLE_OPENMP
    int thread = omp_get_thread_num();
#else
    int thread = 0;
#endif
   currVal[thread] += 1;
  }
  void decrement() {
#ifdef RAJA_ENABLE_OPENMP
    int thread = omp_get_thread_num();
#else
    int thread = 0;
#endif

    currVal[thread] -= 1;
  }

  void reset() {
    for(size_t i = 0; i < numThreads; i++) {
      currVal[0] = low->value();
    }
  }
  
  friend std::ostream& operator<< (std::ostream& s, SymbolicSegment i) {
    s << *(i.low) << " <= " << i.name << " < " << *(i.high);
    return s;
  }

private:
idx_t * currVal;

}; //SymbolicSegment


struct SymbolicBoundAdd : SymbolicBound {
  SymbolicBound *left;
  SymbolicBound * right;

  idx_t value() {
    return left->value() + right->value();
  }

  std::ostream& write(std::ostream & s) {
    s << left << "+" << right;
    return s;
  }
};
struct SymbolicBoundMult : SymbolicBound {
  SymbolicBound * left;
  SymbolicBound * right;

  SymbolicBoundMult(SymbolicBound * l, SymbolicBound * r) : left(l), right(r) {}
  idx_t value() const {
    return left->value() * right->value();
  }

  std::ostream& write(std::ostream & s) {
    s << left << "*" << right;
    return s;
  }

};



SymbolicBound * make_symbolic_bound(idx_t x);
SymbolicBound * make_symbolic_bound(int x);
SymbolicBound * make_symbolic_bound(const SymbolicSegment &x);
SymbolicBound * make_symbolic_bound(SymbolicBound * x);

template <typename T1, typename T2>
SymbolicSegment make_symbolic_segment(const T1 & l, const T2 r, std::string n) {
  
  SymbolicBound * low = make_symbolic_bound(l);
  SymbolicBound * high = make_symbolic_bound(r); 
  return SymbolicSegment(low,high,n);
}

template <typename T1, typename T2>
SymbolicSegment make_symbolic_segment(const T1 & l, const T2 r) {
  static int i = 0;
  i += 1; 
  SymbolicBound * low = make_symbolic_bound(l);
  SymbolicBound * high = make_symbolic_bound(r); 
  std::stringstream name;
  name << "tmpname" << i;
  return SymbolicSegment(low,high,name.str());
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
