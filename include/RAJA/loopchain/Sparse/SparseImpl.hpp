
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
  virtual size_t get_hits() = 0;
  virtual size_t get_misses() = 0;
  virtual void reset_counters() = 0;
  double get_miss_rate() {
    size_t m = get_misses();
    size_t h = get_hits();
    size_t total = m + h;
    return (double) m / (double) total;
  }
  double get_hit_rate() {
    size_t m = get_misses();
    size_t h = get_hits();
    size_t total = m + h;
    return (double) h / (double) total;
  }
  virtual std::vector<IdxType> dim(IdxType) = 0; 
};


}
#endif
