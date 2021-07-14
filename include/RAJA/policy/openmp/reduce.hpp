/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Header file containing RAJA reduction templates for
 *          OpenMP execution.
 *
 *          These methods should work on any platform that supports OpenMP.
 *
 ******************************************************************************
 */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2016-21, Lawrence Livermore National Security, LLC
// and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef RAJA_omp_reduce_HPP
#define RAJA_omp_reduce_HPP

#include "RAJA/config.hpp"

#if defined(RAJA_ENABLE_OPENMP)

#include <memory>
#include <vector>

#include <omp.h>

#include "RAJA/util/types.hpp"

#include "RAJA/pattern/detail/reduce.hpp"
#include "RAJA/pattern/reduce.hpp"

#include "RAJA/policy/openmp/policy.hpp"

namespace RAJA
{

namespace detail
{
template <typename T, typename Reduce>
class ReduceOMP
    : public reduce::detail::BaseCombinable<T, Reduce, ReduceOMP<T, Reduce>>
{
  using Base = reduce::detail::BaseCombinable<T, Reduce, ReduceOMP>;

public:
  using Base::Base;
  //! prohibit compiler-generated default ctor
  ReduceOMP() = delete;

  ~ReduceOMP()
  {
    if (Base::parent) {
#pragma omp critical(ompReduceCritical)
      Reduce()(Base::parent->local(), Base::my_data);
      Base::my_data = Base::identity;
    }
  }


};

template <typename T, typename Reduce>
class ReduceOMPArr
    : public reduce::detail::BaseCombinableArr<T, Reduce, ReduceOMPArr<T, Reduce>>
{
  using Base = reduce::detail::BaseCombinableArr<T, Reduce, ReduceOMPArr>;

public:
  using Base::Base;
  
  //! prohibit compiler-generated default ctor
  ReduceOMPArr() = delete;

  ~ReduceOMPArr()
  {
    //std::cout << "ReduceOmpArr destructor\n";
    if (Base::parent) {
#pragma omp critical(ompReduceCritical)
      Reduce().end_reduce(Base::parent->local(), Base::my_arr, Base::len);
    }
    for(int i = 0; i < Base::len; i++) {
      Base::my_arr[i] = Base::identity;
    }
    free(Base::my_arr);
    Base::len = 0;
    Base::my_arr = NULL;
    
  }

};

template <typename T, typename Reduce>
class ReduceOMPMap
    : public reduce::detail::BaseCombinableMap<T, Reduce, ReduceOMPMap<T, Reduce>>
{
  using Base = reduce::detail::BaseCombinableMap<T, Reduce, ReduceOMPMap>;

public:
  using Base::Base;
  
  //! prohibit compiler-generated default ctor
  ReduceOMPMap() = delete;

  ~ReduceOMPMap()
  {
    //std::cout << "ReduceOmpArr destructor\n";
    if (Base::parent) {
#pragma omp critical(ompReduceCritical)
      //std::cout << "reducing into shared array.\n";
      Reduce().end_reduce(Base::parent->local(), Base::my_map);
    } 
    Base::reset();
  
  }

}; //ReduceOMPMap



}  // namespace detail

RAJA_DECLARE_ALL_REDUCERS(omp_reduce, detail::ReduceOMP)


template <typename T>                                       
  class ReduceSumArr<omp_reduce, T>                                    
      : public reduce::detail::BaseReduceSumArr<T, detail::ReduceOMPArr>    
  {                                                           
  public:                                                     
    using Base = reduce::detail::BaseReduceSumArr<T, detail::ReduceOMPArr>; 
    using Base::Base;                                         
  };

template <typename T>
  class ReduceSumMap<omp_reduce, T>
    : public reduce::detail::BaseReduceSumMap<T, detail::ReduceOMPMap>
{
  public:
    using Base = reduce::detail::BaseReduceSumMap<T, detail::ReduceOMPMap>;
    using Base::Base;
};
///////////////////////////////////////////////////////////////////////////////
//
// Old ordered reductions are included below.
//
///////////////////////////////////////////////////////////////////////////////

namespace detail
{
template <typename T, typename Reduce>
class ReduceOMPOrdered
    : public reduce::detail::
          BaseCombinable<T, Reduce, ReduceOMPOrdered<T, Reduce>>
{
  using Base = reduce::detail::BaseCombinable<T, Reduce, ReduceOMPOrdered>;
  std::shared_ptr<std::vector<T>> data;

public:
  ReduceOMPOrdered() { reset(T(), T()); }

  //! constructor requires a default value for the reducer
  explicit ReduceOMPOrdered(T init_val, T identity_)
  {
    reset(init_val, identity_);
  }

  void reset(T init_val, T identity_)
  {
    Base::reset(init_val, identity_);
    data = std::shared_ptr<std::vector<T>>(
        std::make_shared<std::vector<T>>(omp_get_max_threads(), identity_));
  }

  ~ReduceOMPOrdered()
  {
    Reduce{}((*data)[omp_get_thread_num()], Base::my_data);
    Base::my_data = Base::identity;
  }

  T get_combined() const
  {
    if (Base::my_data != Base::identity) {
      Reduce{}((*data)[omp_get_thread_num()], Base::my_data);
      Base::my_data = Base::identity;
    }

    T res = Base::identity;
    for (size_t i = 0; i < data->size(); ++i) {
      Reduce{}(res, (*data)[i]);
    }
    return res;
  }
};

}  // namespace detail

RAJA_DECLARE_ALL_REDUCERS(omp_reduce_ordered, detail::ReduceOMPOrdered)

template <typename ContentType>
struct DenseReducer<omp_reduce, ContentType> { 
void operator()(ContentType * arr1, const ContentType * arr2, const unsigned & size) {
  #pragma omp critical(ompReduceCriticalDense)
  for(unsigned i = 0; i < size; i++) {
    arr1[i] += arr2[i];
  }
}
}; //DenseReducer<omp>

template <typename ContentType>
struct MapReducer<omp_reduce, ContentType> {
  void operator() (const MapReduction<omp_reduce,ContentType> * out, const MapReduction<omp_reduce,ContentType> * in)   {
    
    if (out->denseContent) {
      #pragma omp critical(ompReduceCriticalMap)
      for(const auto &it : in->mapContent) {
        out->denseContent[it.first] += it.second;
      }
    } else {
      #pragma omp critical(ompReduceCriticalMap)
      for (const auto &it : in->mapContent) {
        (*out)(it.first) += it.second;
      }
    }
  }
};

template <typename ContentType>
struct BlockReducer<omp_reduce, ContentType> {
  void operator() (const BlockReduction<omp_reduce,ContentType> *out, 
                   const BlockReduction<omp_reduce,ContentType> * in) 
  {
    #pragma omp critical(ompReduceCriticalBlock)
    {
    int i, j;
    if (out->nblocks == 1 && out->blocks[0] != NULL) {
      for(i = 0; i < in->nblocks; i++) {
        if(in->blocks[i] != NULL) {
          int startIdx = i * BSIZE;
          if(in->orig->useLocks && in->blocks[i] == &(out->blocks[0][startIdx])) {
            omp_unset_lock(&(out->writelocks[i]));
          } else {
 
            BlockReduction<seq_reduce,ContentType>::addBlock(out->blocks[0] + startIdx, in->blocks[i], i, out->totalsize, false);
            free(in->blocks[i]);
          }
        }
      }
    } else {
      for(i = 0; i < in->nblocks; i++) {
        if(in->blocks[i] != NULL) {
          if(out->blocks[i] != NULL) {
            ContentType *curOut = out->blocks[i];
            ContentType *curIn = in->blocks[i];
            ContentType *origBlock = &(in->orig->blocks[0][i *BSIZE]);
            if(in->orig->useLocks && curIn == origBlock) {
              curOut = in->blocks[i];
              curIn = out->blocks[i];
              out->blocks[i] = in->blocks[i];
            }
            if(in->orig->useLocks && curOut == origBlock) {
              BlockReduction<omp_reduce,ContentType>::addBlock(curOut, curIn, i, out->totalsize, false);
            } else {
              BlockReduction<omp_reduce,ContentType>::addBlock(curOut, curIn, i, out->totalsize, true);
            }
            free(curIn);
          } else {
            out->blocks[i] = in->blocks[i];
          }
        }
      }
    }
    } //omp critical
  }

};


}  // namespace RAJA

#endif  // closing endif for RAJA_ENABLE_OPENMP guard

#endif  // closing endif for header file include guard
