#ifndef _RAJA_SparseSegment_HPP
#define _RAJA_SparseSegment_HPP

#include "RAJA/loopchain/SymbolicSegment/SymbolicSegment.hpp"

namespace RAJA {

template <camp::idx_t NumDims>
struct SparseSymbolicSegmentLead;

template <camp::idx_t NumDims>
struct SparseIteratorWrapper {

  SparseSymbolicSegmentLead<NumDims> * seg;
  using iterator = std::vector<idx_t>::const_iterator; 
  iterator iter;

  SparseIteratorWrapper(SparseSymbolicSegmentLead<NumDims> * s, iterator i) : seg(s), iter(i) {}

  
  idx_t operator[] (idx_t i) {
    seg->currVal = i;
    std::cout << "Set the currVal to " << i << "\n";
    return iter[i];
  }

  idx_t operator - (const SparseIteratorWrapper<NumDims> & other) {
    auto i = iter - other.iter;
    return i;
  } 
};


template <camp::idx_t NumDims>
struct SparseSymbolicSegmentLead {
using iterator = SparseIteratorWrapper<NumDims>;
using value_type = Index_type;
  using IndexType = Index_type;

const SparseView<double, NumDims> & view;
idx_t dimIdx;
std::vector<idx_t> dim;
idx_t currVal = 0;

SparseSymbolicSegmentLead(SparseView<double, NumDims> & v, idx_t d) :
  view(v), dimIdx(d) {
  dim = view.impl->dim(dimIdx);
}


iterator begin() {
  return SparseIteratorWrapper<NumDims>(this, dim.begin());
}

iterator end() {
  return SparseIteratorWrapper<NumDims>(this, dim.end());
}

};


template <camp::idx_t NumDims>
struct SparseIteratorWrapperFollow {

  SparseSymbolicSegmentLead<NumDims> * seg;
  using iterator = std::vector<idx_t>::const_iterator; 
  iterator iter;

  SparseIteratorWrapperFollow(SparseSymbolicSegmentLead<NumDims> * s, iterator i) : seg(s), iter(i) {}

  
  idx_t operator[] (idx_t i) {
    std::cout << "follow wrapper for i, currVal: " << i << " " << seg->currVal << std::endl;
    return iter[i + seg->currVal];
  }

  idx_t operator - (const SparseIteratorWrapperFollow<NumDims> & other) {
    auto i = iter - other.iter;
    return i;
  } 
};




template <camp::idx_t NumDims>
struct SparseSymbolicSegmentFollow {
using iterator = SparseIteratorWrapperFollow<NumDims>;
using value_type = Index_type;
using IndexType = Index_type;

const SparseView<double, NumDims> & view;
idx_t dimIdx;
std::vector<idx_t> dim;
SparseSymbolicSegmentLead<NumDims> * leadSeg;

SparseSymbolicSegmentFollow(SparseView<double, NumDims> & v, idx_t d, SparseSymbolicSegmentLead<NumDims> * l) :
  view(v), dimIdx(d), leadSeg(l) {
  dim = view.impl->dim(dimIdx);
}


iterator begin() const {
  std::cout << "begin leadSeg.currVal:" << leadSeg->currVal << "\n";
  return SparseIteratorWrapperFollow<NumDims>(leadSeg, dim.begin());
}

iterator end() const {
  std::cout << "end leadSeg.currVal:" << leadSeg->currVal << "\n";
  return SparseIteratorWrapperFollow<NumDims>(leadSeg, dim.begin() + 1);
}

};



} //namespace RAJA


template <camp::idx_t NumDims>
struct std::iterator_traits<RAJA::SparseIteratorWrapper<NumDims>> {
  typedef RAJA::Index_type value_type;
  typedef RAJA::Index_type difference_type;
  typedef std::random_access_iterator_tag iterator_category;
  typedef RAJA::Index_type* pointer;
  typedef RAJA::Index_type& reference;
};
template <camp::idx_t NumDims>
struct std::iterator_traits<RAJA::SparseIteratorWrapperFollow<NumDims>> {
  typedef RAJA::Index_type value_type;
  typedef RAJA::Index_type difference_type;
  typedef std::random_access_iterator_tag iterator_category;
  typedef RAJA::Index_type* pointer;
  typedef RAJA::Index_type& reference;
};



#endif
