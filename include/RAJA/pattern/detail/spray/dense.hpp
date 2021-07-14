#ifndef RAJA_DETAIL_REDUCE_SPRAY_DENSE_HPP
#define RAJA_DETAIL_REDUCE_SPRAY_DENSE_HPP

namespace RAJA {



template <typename ReductionPolicy, typename ContentType>
class DenseReducer {};


template <typename ReductionPolicy, typename ContentType> 
class DenseReduction 
{
public:
  ContentType * content = nullptr;
  DenseReduction<ReductionPolicy, ContentType> const * parent = nullptr;
  

  unsigned size = 0;


  DenseReduction() {}

  DenseReduction(unsigned size_, ContentType *orig) 
    : content(orig), size(size_) 
  {
  }

  DenseReduction(const DenseReduction &other) 
    : parent(other.parent ? other.parent : &other),
      size(other.size)
  {
    content = new ContentType[size];
    for(unsigned i = 0; i < size; i++) {content[i] = 0;}
  }

  ~DenseReduction() 
  {
    if(parent && content) {
      DenseReducer<ReductionPolicy,ContentType>{}(parent->content, content, size);
      size = 0;
      free(content);
      content = nullptr;
    }
  }

  ContentType &operator()(int idx) const {
    return content[idx];
  }


}; // DenseReduction

} // namespace RAJA





#endif
