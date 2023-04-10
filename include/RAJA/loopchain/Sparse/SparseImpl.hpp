
#ifndef _RAJA_SPARSEIMPL_HPP
#define _RAJA_SPARSEIMPL_HPP

#include <cassert>

namespace RAJA {

template <typename ElmType,size_t NumDims>
struct SparseImpl {

using IdxType = camp::idx_t;
virtual ElmType & operator() (std::vector<IdxType> v) = 0;
  
  virtual std::vector<IdxType> get_permutation() = 0;  
  virtual void set_permutation(std::vector<IdxType>) = 0;  
  virtual void* get_data() = 0; 

  virtual std::vector<IdxType> dim(IdxType) = 0; 
};


}
#endif
