#ifndef _RAJA_LCENUMERATE_HPP
#define _RAJA_LCENUMERATE_HPP

#include "RAJA/loopchain/SymExec.hpp"
#include <algorithm>
namespace RAJA {


template <typename ViewType>
auto enumerate_perms(ViewType v) {

  constexpr camp::idx_t NumDims = v.get_layout().n_dims;

  using ArrayType = std::array<camp::idx_t, NumDims>;

  std::vector<ArrayType> perms;

  ArrayType currPerm;

   
  for(int i = 0; i < NumDims; i++) {
    currPerm[i] = layout_to_perm(v.get_layout())[i];
  }    

    std::sort(currPerm.begin(), currPerm.end());

  do {
    ArrayType perm;
    for(int i = 0; i < NumDims; i++) {
      perm[i] = currPerm[i];
    }
    perms.push_back(perm);
  } while (std::next_permutation(currPerm.begin(), currPerm.end()));
  
  
  return perms;  
  
}

template <camp::idx_t NumDims>
auto enumerate_perms() {

  using ArrayType = std::array<camp::idx_t, NumDims>;

  std::vector<ArrayType> perms;

  ArrayType currPerm;

   
  for(int i = 0; i < NumDims; i++) {
    currPerm[i] = i;
  }    

    std::sort(currPerm.begin(), currPerm.end());

  do {
    ArrayType perm;
    for(int i = 0; i < NumDims; i++) {
      perm[i] = currPerm[i];
    }
    perms.push_back(perm);
  } while (std::next_permutation(currPerm.begin(), currPerm.end()));
  
  
  return perms;  

}





} //namespace RAJA

#endif