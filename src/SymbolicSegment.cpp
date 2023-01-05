#include "RAJA/loopchain/SymbolicSegment/SymbolicSegment.hpp"

namespace RAJA {

idx_t SymbolicBoundSegment::value() const {
    return seg->value();
}
std::ostream & SymbolicBoundSegment::write(std::ostream & s) {
  s << seg->name;
  return s;
}

SymbolicBound * SymbolicSegment::operator * (idx_t c) {
  return new SymbolicBoundMult(new SymbolicBoundSegment(*this), new SymbolicBoundConstant(c));
}
SymbolicBound * SymbolicSegment::operator + (idx_t c) {
  return new SymbolicBoundAdd(new SymbolicBoundSegment(*this), new SymbolicBoundConstant(c));
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

}
