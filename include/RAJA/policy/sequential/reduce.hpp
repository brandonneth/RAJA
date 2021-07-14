/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Header file containing RAJA reduction templates for
 *          sequential execution.
 *
 *          These methods should work on any platform.
 *
 ******************************************************************************
 */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2016-21, Lawrence Livermore National Security, LLC
// and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef RAJA_sequential_reduce_HPP
#define RAJA_sequential_reduce_HPP

#include "RAJA/config.hpp"

#include "RAJA/internal/MemUtils_CPU.hpp"

#include "RAJA/pattern/detail/reduce.hpp"
#include "RAJA/pattern/reduce.hpp"

#include "RAJA/policy/sequential/policy.hpp"

#include "RAJA/util/types.hpp"

namespace RAJA
{

namespace detail
{
template <typename T, typename Reduce>
class ReduceSeq
    : public reduce::detail::BaseCombinable<T, Reduce, ReduceSeq<T, Reduce>>
{
  using Base = reduce::detail::BaseCombinable<T, Reduce, ReduceSeq<T, Reduce>>;

public:
  //! prohibit compiler-generated default ctor
  ReduceSeq() = delete;

  using Base::Base;
};
template <typename T, typename Reduce>
class ReduceSeqArr
    : public reduce::detail::BaseCombinableArr<T, Reduce, ReduceSeqArr<T, Reduce>>
{
  using Base = reduce::detail::BaseCombinableArr<T, Reduce, ReduceSeqArr<T, Reduce>>;

public:
  //! prohibit compiler-generated default ctor
  ReduceSeqArr() = delete;

  using Base::Base;
};

}  // namespace detail

RAJA_DECLARE_ALL_REDUCERS(seq_reduce, detail::ReduceSeq)

//RAJA_DECLARE_REDUCER(SumArr, seq_reduce, detail::ReduceSeq)

template <typename T>
class ReduceSumArr<seq_reduce, T> : public reduce::detail::BaseReduceSumArr<T, detail::ReduceSeqArr>
{
public:
  using Base = reduce::detail::BaseReduceSumArr<T, detail::ReduceSeqArr>;
  using Base::Base;
};


template <typename ContentType>
class DenseReducer<seq_reduce, ContentType> {

public:
  void operator() (ContentType * arr1, ContentType * arr2, const unsigned & size) {
    std::cout << "Sequential dense reducer\n";

    for(int i = 0; i < size; i++) {
      arr1[i] += arr2[i];
    }
  }

};

template <typename ContentType>
struct MapReducer<seq_reduce, ContentType> {
  void operator() (const MapReduction<seq_reduce,ContentType> * out, const MapReduction<seq_reduce,ContentType> * in)   {

    if (out->denseContent) {
      for(const auto &it : in->mapContent) {
        out->denseContent[it.first] += it.second;
      }
    } else {
      for (const auto &it : in->mapContent) {
        (*out)(it.first) += it.second;
      }
    }
  }
};

template <typename ContentType>
struct BlockReducer<seq_reduce, ContentType> {
  void operator() (const BlockReduction<seq_reduce,ContentType> *out, 
                   const BlockReduction<seq_reduce,ContentType> * in) 
  {
    int i, j;
    if (out->nblocks == 1 && out->blocks[0] != NULL) {
      for(i = 0; i < in->nblocks; i++) {
        if(in->blocks[i] != NULL) {
          int startIdx = i * BSIZE;
          BlockReduction<seq_reduce,ContentType>::addBlock(out->blocks[0] + startIdx, in->blocks[i], i, out->totalsize, false);
          free(in->blocks[i]);
        }
      }
    } else {
      for(i = 0; i < in->nblocks; i++) {
        if(in->blocks[i] != NULL) {
          if(out->blocks[i] != NULL) {
            ContentType *curOut = out->blocks[i];
            ContentType *curIn = in->blocks[i];
            ContentType *origBlock = &(in->orig->blocks[0][i *BSIZE]);
            BlockReduction<seq_reduce,ContentType>::addBlock(curOut, curIn, i, out->totalsize, true);
            free(curIn);
          } else {
            out->blocks[i] = in->blocks[i];
          }
        }
      }
    }
  }
};


}  // namespace RAJA

#endif  // closing endif for header file include guard
