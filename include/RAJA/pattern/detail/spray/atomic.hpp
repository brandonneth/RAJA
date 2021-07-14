#ifndef RAJA_REDUCE_SPRAY_ATOMIC_HPP
#define RAJA_REDUCE_SPRAY_ATOMIC_HPP

namespace RAJA {



template <typename ContentType> class AtomicScalar {
public:
  AtomicScalar(ContentType *val) { this->val = val; }

  void operator+=(ContentType rhs) {
#pragma omp atomic update
    *(this->val) += rhs;
  }

private:
  ContentType *val;
};


template <typename ReductionPolicy, typename ContentType> 
class AtomicReduction 
{
public:
  ContentType * content = nullptr;
  AtomicReduction<ReductionPolicy, ContentType> const * parent = nullptr;
  

  unsigned size = 0;


  AtomicReduction() {}

  AtomicReduction(unsigned size_, ContentType *orig) 
    : content(orig), size(size_) 
  {
  }

  AtomicReduction(const AtomicReduction &other) 
    : content(other.content),
      size(other.size)
  {
  }

  ~AtomicReduction() 
  {
  }

  AtomicScalar<ContentType> operator()(int idx) const {
    return AtomicScalar<ContentType>(content + idx);
  }


}; // DenseReduction

} // namespace RAJA





#endif
