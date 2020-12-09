/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Header file for tile wrapper and iterator.
 *
 ******************************************************************************
 */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Copyright (c) 2016-20, Lawrence Livermore National Security, LLC
// and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifndef RAJA_pattern_kernel_OverlappedTile_HPP
#define RAJA_pattern_kernel_OverlappedTile_HPP

#include "RAJA/config.hpp"

#include <iostream>
#include <type_traits>

#include "camp/camp.hpp"
#include "camp/concepts.hpp"
#include "camp/tuple.hpp"

#include "RAJA/pattern/kernel/Tile.hpp"
#include "RAJA/pattern/kernel/internal.hpp"
#include "RAJA/util/macros.hpp"
#include "RAJA/util/types.hpp"

namespace RAJA
{
namespace statement
{


/*!
 * A RAJA::kernel statement that implements a tiling (or blocking) loop.
 *
 */
template <camp::idx_t ArgumentId,
          typename TilePolicy,
          typename ExecPolicy,
          typename... EnclosedStmts>
struct OverlappedTile : public internal::Statement<ExecPolicy, EnclosedStmts...> {
  using tile_policy_t = TilePolicy;
  using exec_policy_t = ExecPolicy;
};


}  // end namespace statement

namespace internal
{

/*!
 * A generic RAJA::kernel forall_impl tile wrapper for statement::For
 * Assigns the tile segment to segment ArgumentId
 *
 */
template <camp::idx_t ArgumentId, typename Data, typename Types, typename... EnclosedStmts>
struct OverlappedTileWrapper : public GenericWrapper<Data, Types, EnclosedStmts...> {

  using Base = GenericWrapper<Data, Types, EnclosedStmts...>;
  using Base::Base;
  using privatizer = NestedPrivatizer<OverlappedTileWrapper>;

  template <typename InSegmentIndexType>
  RAJA_INLINE void operator()(InSegmentIndexType si)
  {
    // Assign the tile's segment to the tuple
    camp::get<ArgumentId>(Base::data.segment_tuple) = si.s;
    //Base::data.template assign_offset<ArgumentId>(si);
    // Execute enclosed statements
    //std::cout << "calling Base::exec()\n";
    Base::exec();
  }
};


template <typename Iterable>
struct IterableOverlappedTiler {
  using value_type = camp::decay<Iterable>;

  struct iterate
  {
    value_type s;
    Index_type i;
  };

  class iterator
  {
    // NOTE: this must be held by value for NVCC support, *even on the host*
    const IterableOverlappedTiler itiler;
    const Index_type block_id;
  public:
    using value_type = iterate;
    using difference_type = camp::idx_t;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::random_access_iterator_tag;

    RAJA_HOST_DEVICE
    RAJA_INLINE
    constexpr iterator(IterableOverlappedTiler const &itiler_, Index_type block_id_)
        : itiler{itiler_}, block_id{block_id_}    {
    }

    RAJA_HOST_DEVICE
    RAJA_INLINE
    value_type operator*()
    {
      auto start = block_id * itiler.block_size;
      return iterate{itiler.it.slice(start, itiler.block_size), block_id};
    }

    RAJA_HOST_DEVICE
    RAJA_INLINE difference_type operator-(const iterator &rhs) const
    {
      return static_cast<difference_type>(block_id) -
             static_cast<difference_type>(rhs.block_id);
    }

    RAJA_HOST_DEVICE
    RAJA_INLINE iterator operator-(const difference_type &rhs) const
    {
      return iterator(itiler, block_id - rhs);
    }

    RAJA_HOST_DEVICE
    RAJA_INLINE iterator operator+(const difference_type &rhs) const
    {
      return iterator(itiler,
                      block_id + rhs >= itiler.num_blocks ? itiler.num_blocks
                                                          : block_id + rhs);
    }

    RAJA_HOST_DEVICE
    RAJA_INLINE value_type operator[](difference_type rhs) const
    {
      return *((*this) + rhs);
    }

    RAJA_HOST_DEVICE
    RAJA_INLINE bool operator!=(const IterableOverlappedTiler &rhs) const
    {
      return block_id != rhs.block_id;
    }

    RAJA_HOST_DEVICE
    RAJA_INLINE bool operator<(const IterableOverlappedTiler &rhs) const
    {
      return block_id < rhs.block_id;
    }
  };

  RAJA_HOST_DEVICE
  RAJA_INLINE
  IterableOverlappedTiler(const Iterable &it_, camp::idx_t block_size_)
      : it{it_}, block_size{block_size_}  {
    using std::begin;
    using std::distance;
    using std::end;
    dist = it.end() - it.begin();  // distance(begin(it), end(it));
    num_blocks = dist / block_size;
    // if (dist % block_size) num_blocks += 1;
    if (dist - num_blocks * block_size > 0) {
      num_blocks += 1;
    }
  }

  RAJA_HOST_DEVICE
  RAJA_INLINE
  iterator begin() const { return iterator(*this, 0); }

  RAJA_HOST_DEVICE
  RAJA_INLINE
  iterator end() const { return iterator(*this, num_blocks); }

  value_type it;
  camp::idx_t block_size;
  camp::idx_t num_blocks;
  camp::idx_t dist;
};

/*!
 * A generic RAJA::kernel forall_impl executor for statement::Tile
 *
 *
 */
template <camp::idx_t ArgumentId,
          typename TPol,
          typename EPol,
          typename... EnclosedStmts,
          typename Types>
struct StatementExecutor<
    statement::OverlappedTile<ArgumentId, TPol, EPol, EnclosedStmts...>, Types> {


  template <typename Data>
  static RAJA_INLINE void exec(Data &data)
  {
    //std::cout << "StatementExecutor<OverlappedTile>\n";
    //std::cout << "Executing statement executor for overlappedtile\n";
    // Get the segment we are going to tile
    auto const &segment = camp::get<ArgumentId>(data.segment_tuple);

    // Get the tiling policies chunk size
    auto chunk_size = TPol::chunk_size;
    
    auto tileSizes = data.tileSizes;
    auto tileSize = chunk_size; 
    if (tileSizes.size() > ArgumentId) {
      //std::cout << "tileSize: " << tileSizes.at(ArgumentId) << "\n";
      tileSize = tileSizes.at(ArgumentId);
    } else {
      tileSize = chunk_size;
    }
    // Create a tile iterator, needs to survive until the forall is
    // done executing.
    IterableOverlappedTiler<decltype(segment)> tiled_iterable(segment, tileSize);

    // Wrap in case forall_impl needs to thread_privatize
    OverlappedTileWrapper<ArgumentId, Data, Types,
                EnclosedStmts...> tile_wrapper(data);

    auto r = resources::get_resource<EPol>::type::get_default();
    // Loop over tiles, executing enclosed statement list
    forall_impl(r, EPol{}, tiled_iterable, tile_wrapper);

    // Set range back to original values
    camp::get<ArgumentId>(data.segment_tuple) = tiled_iterable.it;
  }
};

}  // end namespace internal
}  // end namespace RAJA

#endif /* RAJA_pattern_kernel_HPP */