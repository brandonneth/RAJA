/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   RAJA header file defining Layout, a N-dimensional index calculator
 *          with offset indices
 *
 ******************************************************************************
 */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2016-23, Lawrence Livermore National Security, LLC
// and RAJA project contributors. See the RAJA/LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef RAJA_MORTONLAYOUT_HPP
#define RAJA_MORTONLAYOUT_HPP

#include "RAJA/config.hpp"

#include <array>
#include <limits>

#include "camp/camp.hpp"

#include "RAJA/index/IndexValue.hpp"

#include "RAJA/util/Permutations.hpp"
#include "RAJA/util/PermutedLayout.hpp"

namespace RAJA
{

namespace internal
{

template <typename Range, typename IdxLin>
struct MortonLayout_impl;

template <camp::idx_t... RangeInts, typename IdxLin>
struct MortonLayout_impl<camp::idx_seq<RangeInts...>, IdxLin> {
  using Self = MortonLayout_impl<camp::idx_seq<RangeInts...>, IdxLin>;
  using IndexRange = camp::idx_seq<RangeInts...>;
  using IndexLinear = IdxLin;
  using Base = RAJA::detail::LayoutBase_impl<IndexRange, IdxLin>;
  Base base_;

  static constexpr camp::idx_t stride_one_dim = Base::stride_one_dim;

  static constexpr size_t n_dims = sizeof...(RangeInts);
  IdxLin dimSizes[n_dims]={0}; //If not specified set to zero
  IdxLin tileSizes[n_dims]={0}; //If not specified set to zero

  constexpr RAJA_INLINE MortonLayout_impl(
      std::array<IdxLin, sizeof...(RangeInts)> _dimSizes,
      std::array<IdxLin, sizeof...(RangeInts)> _tileSizes)
      : base_{_tileSizes[RangeInts]...},
        dimSizes{_dimSizes[RangeInts]...},
        tileSizes{_tileSizes[RangeInts]...}
  {
  }

  constexpr RAJA_INLINE RAJA_HOST_DEVICE MortonLayout_impl(Self const& c)
      : base_(c.base_), dimSizes{c.dimSizes[RangeInts]...}
  {
  }
  template<camp::idx_t N, typename Idx>
  RAJA_INLINE RAJA_HOST_DEVICE void BoundsCheckError(Idx idx) const
  {
    printf("Error at index %d, value %ld is not within bounds [%ld, %ld] \n",
           static_cast<int>(N), static_cast<long int>(idx),
           static_cast<long int>(dimSizes[N]), static_cast<long int>(dimSizes[N] + base_.sizes[N] - 1));
    RAJA_ABORT_OR_THROW("Out of bounds error \n");
  }

  template <camp::idx_t N>
  RAJA_INLINE RAJA_HOST_DEVICE void BoundsCheck() const
  {
  }

  template <camp::idx_t N, typename Idx, typename... Indices>
  RAJA_INLINE RAJA_HOST_DEVICE void BoundsCheck(Idx idx, Indices... indices) const
  {
    if(!(dimSizes[N] <=idx && idx < dimSizes[N] + base_.sizes[N]))
    {
      BoundsCheckError<N>(idx);
    }
    RAJA_UNUSED_VAR(idx);
    BoundsCheck<N+1>(indices...);
  }

  unsigned int interleave(unsigned int x, unsigned int y) const {
    
    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;

    y = (y | (y << 8)) & 0x00FF00FF;
    y = (y | (y << 4)) & 0x0F0F0F0F;
    y = (y | (y << 2)) & 0x33333333;
    y = (y | (y << 1)) & 0x55555555;

    unsigned int z = (x << 1) | y;
    return z;
  }

  template <typename... Indices>
  RAJA_INLINE RAJA_HOST_DEVICE RAJA_BOUNDS_CHECK_constexpr IdxLin operator()(
      Indices... indices) const
  {
#if defined (RAJA_BOUNDS_CHECK_INTERNAL)
    BoundsCheck<0>(indices...);
#endif
    auto tileIndex = interleave((indices / tileSizes[RangeInts])...);
    auto tileStartOffset = tileIndex * (tileSizes[RangeInts] * ...);

    auto inTileOffset = base_((indices % tileSizes[RangeInts])...);
   
    return tileStartOffset + inTileOffset;
  }

  template <typename... Indices>
  RAJA_INLINE RAJA_HOST_DEVICE void toIndices(IdxLin linear_index,
                                              Indices &&... indices) const
  {
    
  }


  constexpr RAJA_INLINE RAJA_HOST_DEVICE
  MortonLayout_impl(const Layout<sizeof...(RangeInts), IdxLin>& rhs)
      : base_{rhs}
  {
  }

  RAJA_INLINE RAJA_HOST_DEVICE constexpr IdxLin size() const
  {
    return base_.size();
  }

  RAJA_INLINE RAJA_HOST_DEVICE constexpr IdxLin size_noproj() const
  {
    return base_.size_noproj();
  }

  template<camp::idx_t DIM>
  RAJA_INLINE
  RAJA_HOST_DEVICE
  constexpr
  IndexLinear get_dim_stride() const {
    return base_.get_dim_stride();
  }

  template<camp::idx_t DIM>
  RAJA_INLINE
  RAJA_HOST_DEVICE
  constexpr
  IndexLinear get_dim_size() const {
    return base_.get_dim_size();
  }

};

}  // namespace internal

template <size_t n_dims = 1, typename IdxLin = Index_type>
struct MortonLayout
    : public internal::MortonLayout_impl<camp::make_idx_seq_t<n_dims>, IdxLin> {
  using Base =
      internal::MortonLayout_impl<camp::make_idx_seq_t<n_dims>, IdxLin>;

  using internal::MortonLayout_impl<camp::make_idx_seq_t<n_dims>,
                                    IdxLin>::MortonLayout_impl;

  constexpr RAJA_INLINE RAJA_HOST_DEVICE MortonLayout(
      const internal::MortonLayout_impl<camp::make_idx_seq_t<n_dims>, IdxLin>&
          rhs)
      : Base{rhs}
  {
  }
};


template <size_t n_dims, typename IdxLin = Index_type>
auto make_morton_layout(std::array<IdxLin,n_dims> dimSizes, std::array<IdxLin,n_dims> tileSizes)
{
  
  return MortonLayout<n_dims, IdxLin>{dimSizes, tileSizes};
}


}  // namespace RAJA

#endif
