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
    : public reduce::detail::BaseCombinable<T, Reduce, ReduceSeqArr<T, Reduce>>
{
  using Base = reduce::detail::BaseCombinable<T, Reduce, ReduceSeqArr<T, Reduce>>;

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
}  // namespace RAJA

#endif  // closing endif for header file include guard
