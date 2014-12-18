/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Header file containing RAJA index set iteration template methods 
 *          using for Intel Cilk Plus execution.
 *
 *          These methods work only on platforms that support Cilk Plus. 
 *
 * \author  Rich Hornung, Center for Applied Scientific Computing, LLNL
 * \author  Jeff Keasler, Applications, Simulations And Quality, LLNL 
 *
 ******************************************************************************
 */

#ifndef RAJA_forall_cilk_icc_HXX
#define RAJA_forall_cilk_icc_HXX

#include "config.hxx"

#include "int_datatypes.hxx"

#include "RAJAVec.hxx"

#include "execpolicy.hxx"

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>


namespace RAJA {

//
//////////////////////////////////////////////////////////////////////
//
// Function templates that iterate over range index sets.
//
//////////////////////////////////////////////////////////////////////
//

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over index range.
 *
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall(cilk_for_exec,
            const Index_type begin, const Index_type end, 
            LOOP_BODY loop_body)
{
   RAJA_FT_BEGIN ;

   cilk_for ( Index_type ii = begin ; ii < end ; ++ii ) {
      loop_body( ii );
   }

   RAJA_FT_END ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over index range, including index count.
 *
 *         NOTE: lambda loop body requires two args (icount, index).
 *
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall_Icount(cilk_for_exec,
                   const Index_type begin, const Index_type end,
                   const Index_type icount,
                   LOOP_BODY loop_body)
{
   const Index_type loop_end = end - begin + 1;

   RAJA_FT_BEGIN ;

   cilk_for ( Index_type ii = 0 ; ii < loop_end ; ++ii ) {
      loop_body( ii+icount, ii+begin );
   }

   RAJA_FT_END ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  iteration over index range set object.
 *
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall(cilk_for_exec,
            const RangeSegment& is,
            LOOP_BODY loop_body)
{
   const Index_type begin = is.getBegin();
   const Index_type end   = is.getEnd();

   RAJA_FT_BEGIN ;

   cilk_for ( Index_type ii = begin ; ii < end ; ++ii ) {
      loop_body( ii );
   }

   RAJA_FT_END ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over index range set object, 
 *         including index count.
 *
 *         NOTE: lambda loop body requires two args (icount, index).
 *
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall_Icount(cilk_for_exec,
                   const RangeSegment& is,
                   const Index_type icount,
                   LOOP_BODY loop_body)
{
   const Index_type begin = is.getBegin();
   const Index_type loop_end = is.getEnd() - begin + 1;

   RAJA_FT_BEGIN ;

   cilk_for ( Index_type ii = 0 ; ii < loop_end ; ++ii ) {
      loop_body( ii+icount, ii+begin );
   }

   RAJA_FT_END ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  minloc reduction over index range.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_minloc(cilk_for_exec,
                   const Index_type begin, const Index_type end,
                   T* min, Index_type *loc,
                   LOOP_BODY loop_body)
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> min_tmp(nworkers);
   RAJAVec<Index_type> loc_tmp(nworkers);

   RAJA_FT_BEGIN ;

   for ( int i = 0; i < nworkers; ++i ) {
       min_tmp[i] = *min ;
       loc_tmp[i] = *loc ;
   }

   cilk_for ( Index_type ii = begin ; ii < end ; ++ii ) {
      loop_body( ii, &min_tmp[__cilkrts_get_worker_number()], 
                     &loc_tmp[__cilkrts_get_worker_number()] );
   }

   for ( int i = 1; i < nworkers; ++i ) {
      if ( min_tmp[i] < min_tmp[0] ) {
         min_tmp[0] = min_tmp[i];
         loc_tmp[0] = loc_tmp[i];
      }
   }

   RAJA_FT_END ;

   *min = min_tmp[0] ;
   *loc = loc_tmp[0] ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  minloc reduction over range index set object.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_minloc(cilk_for_exec,
                   const RangeSegment& is,
                   T* min, Index_type *loc,
                   LOOP_BODY loop_body)
{
   forall_minloc(cilk_for_exec(),
                 is.getBegin(), is.getEnd(),
                 min, loc,
                 loop_body);
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  maxloc reduction over index range.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_maxloc(cilk_for_exec,
                   const Index_type begin, const Index_type end,
                   T* max, Index_type *loc,
                   LOOP_BODY loop_body)
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> max_tmp(nworkers);
   RAJAVec<Index_type> loc_tmp(nworkers);

   RAJA_FT_BEGIN ;

   for ( int i = 0; i < nworkers; ++i ) {
       max_tmp[i] = *max ;
       loc_tmp[i] = *loc ;
   }

   cilk_for ( Index_type ii = begin ; ii < end ; ++ii ) {
      loop_body( ii, &max_tmp[__cilkrts_get_worker_number()],
                     &loc_tmp[__cilkrts_get_worker_number()] );
   }  

   for ( int i = 1; i < nworkers; ++i ) {
      if ( max_tmp[i] > max_tmp[0] ) {
         max_tmp[0] = max_tmp[i];
         loc_tmp[0] = loc_tmp[i];
      }
   }

   RAJA_FT_END ;

   *max = max_tmp[0] ;
   *loc = loc_tmp[0] ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  maxloc reduction over range index set object.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_maxloc(cilk_for_exec,
                   const RangeSegment& is,
                   T* max, Index_type *loc,
                   LOOP_BODY loop_body)
{
   forall_maxloc(cilk_for_exec(),
                 is.getBegin(), is.getEnd(),
                 max, loc,
                 loop_body);
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  sum reduction over index range.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_sum(cilk_for_exec,
                const Index_type begin, const Index_type end,
                T* sum,
                LOOP_BODY loop_body)
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> sum_tmp(nworkers);

   RAJA_FT_BEGIN ;

   for ( int i = 0; i < nworkers; ++i ) {
       sum_tmp[i] = 0 ;
   }

   cilk_for ( Index_type ii = begin ; ii < end ; ++ii ) {
      loop_body( ii, &sum_tmp[__cilkrts_get_worker_number()] );
   }

   RAJA_FT_END ;

   for ( int i = 0; i < nworkers; ++i ) {
      *sum += sum_tmp[i];
   }
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  sum reduction over range index set object.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_sum(cilk_for_exec,
                const RangeSegment& is,
                T* sum,
                LOOP_BODY loop_body)
{
   forall_sum(cilk_for_exec(),
              is.getBegin(), is.getEnd(),
              sum,
              loop_body);
}


//
//////////////////////////////////////////////////////////////////////
//
// Function templates that iterate over range index sets with stride.
//
//////////////////////////////////////////////////////////////////////
//

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over index range with stride.
 *         
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall(cilk_for_exec,
            const Index_type begin, const Index_type end, 
            const Index_type stride,
            LOOP_BODY loop_body)
{                    
   RAJA_FT_BEGIN ;

   cilk_for ( Index_type ii = begin ; ii < end ; ii += stride ) {
      loop_body( ii );
   }

   RAJA_FT_END ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over index range with stride,
 *         including index count.
 *
 *         NOTE: lambda loop body requires two args (icount, index).
 *
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall_Icount(cilk_for_exec,
                   const Index_type begin, const Index_type end,
                   const Index_type stride,
                   const Index_type icount,
                   LOOP_BODY loop_body)
{
   const Index_type loop_end = (end-begin)/stride + 1;

   RAJA_FT_BEGIN ;

   cilk_for ( Index_type ii = 0 ; ii < loop_end ; ++ii ) {
      loop_body( ii+icount, begin + ii*stride );
   }

   RAJA_FT_END ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over range index set with stride object.
 *
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall(cilk_for_exec,
            const RangeStrideSegment& is,
            LOOP_BODY loop_body)
{
   const Index_type begin  = is.getBegin();
   const Index_type end    = is.getEnd();
   const Index_type stride = is.getStride();

   RAJA_FT_BEGIN ;

   cilk_for ( Index_type ii = begin ; ii < end ; ii += stride ) {
      loop_body( ii );
   }

   RAJA_FT_END ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over range index set with stride object,
 *         including index count.
 *
 *         NOTE: lambda loop body requires two args (icount, index).
 *
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall_Icount(cilk_for_exec,
                   const RangeStrideSegment& is,
                   const Index_type icount,
                   LOOP_BODY loop_body)
{
   const Index_type begin    = is.getBegin();
   const Index_type stride   = is.getStride();
   const Index_type loop_end = (is.getEnd()-begin)/stride + 1;

   RAJA_FT_BEGIN ;

   cilk_for ( Index_type ii = 0 ; ii < loop_end ; ++ii ) {
      loop_body( ii+icount, begin + ii*stride );
   }

   RAJA_FT_END ;
}


/*!
 ******************************************************************************
 *
 * \brief  cilk_for  minloc reduction over index range with stride.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_minloc(cilk_for_exec,
                   const Index_type begin, const Index_type end, 
                   const Index_type stride,
                   T* min, Index_type *loc,
                   LOOP_BODY loop_body)
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> min_tmp(nworkers);
   RAJAVec<Index_type> loc_tmp(nworkers);

   RAJA_FT_BEGIN ;

   for ( int i = 0; i < nworkers; ++i ) {
       min_tmp[i] = *min ;
       loc_tmp[i] = *loc ;
   }

   cilk_for ( Index_type ii = begin ; ii < end ; ii += stride ) {
      loop_body( ii, &min_tmp[__cilkrts_get_worker_number()],
                     &loc_tmp[__cilkrts_get_worker_number()] );
   }

   for ( int i = 1; i < nworkers; ++i ) {
      if ( min_tmp[i] < min_tmp[0] ) {
         min_tmp[0] = min_tmp[i];
         loc_tmp[0] = loc_tmp[i];
      }
   }

   RAJA_FT_END ;

   *min = min_tmp[0] ;
   *loc = loc_tmp[0] ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for minloc reduction over range index set with stride object.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_minloc(cilk_for_exec,
                   const RangeStrideSegment& is,
                   T* min, Index_type *loc,
                   LOOP_BODY loop_body)
{
   forall_minloc(cilk_for_exec(),
                 is.getBegin(), is.getEnd(), is.getStride(),
                 min, loc,
                 loop_body);
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  maxloc reduction over index range with stride.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_maxloc(cilk_for_exec,
                   const Index_type begin, const Index_type end,
                   const Index_type stride,
                   T* max, Index_type *loc,
                   LOOP_BODY loop_body)
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> max_tmp(nworkers);
   RAJAVec<Index_type> loc_tmp(nworkers);

   RAJA_FT_BEGIN ;

   for ( int i = 0; i < nworkers; ++i ) {
       max_tmp[i] = *max ;
       loc_tmp[i] = *loc ;
   }

   cilk_for ( Index_type ii = begin ; ii < end ; ii += stride ) {
      loop_body( ii, &max_tmp[__cilkrts_get_worker_number()],
                     &loc_tmp[__cilkrts_get_worker_number()] );
   }

   for ( int i = 1; i < nworkers; ++i ) {
      if ( max_tmp[i] > max_tmp[0] ) {
         max_tmp[0] = max_tmp[i];
         loc_tmp[0] = loc_tmp[i];
      }
   }

   RAJA_FT_END ;

   *max = max_tmp[0] ;
   *loc = loc_tmp[0] ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for maxloc reduction over range index set with stride object.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_maxloc(cilk_for_exec,
                   const RangeStrideSegment& is,
                   T* max, Index_type *loc,
                   LOOP_BODY loop_body)
{
   forall_maxloc(cilk_for_exec(),
                 is.getBegin(), is.getEnd(), is.getStride(),
                 max, loc,
                 loop_body);
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  sum reduction over index range with stride.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_sum(cilk_for_exec,
                const Index_type begin, const Index_type end,
                const Index_type stride,
                T* sum,
                LOOP_BODY loop_body)
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> sum_tmp(nworkers);

   RAJA_FT_BEGIN ;

   for ( int i = 0; i < nworkers; ++i ) {
      sum_tmp[i] = 0 ;
   }

   cilk_for ( Index_type ii = begin ; ii < end ; ii += stride ) {
      loop_body( ii, &sum_tmp[__cilkrts_get_worker_number()] );
   }

   RAJA_FT_END ;

   for ( int i = 0; i < nworkers; ++i ) {
      *sum += sum_tmp[i];
   }
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  sum reduction over range index set with stride object.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_sum(cilk_for_exec,
                const RangeStrideSegment& is,
                T* sum,
                LOOP_BODY loop_body)
{
   forall_sum(cilk_for_exec(),
              is.getBegin(), is.getEnd(), is.getStride(),
              sum,
              loop_body);
}



//
//////////////////////////////////////////////////////////////////////
//
// Function templates that iterate over unstructured index sets.
//
//////////////////////////////////////////////////////////////////////
//

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over indices in indirection array.
 *
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall(cilk_for_exec,
            const Index_type* __restrict__ idx, const Index_type len,
            LOOP_BODY loop_body)
{
   RAJA_FT_BEGIN ;

   cilk_for ( Index_type k = 0 ; k < len ; ++k ) {
      loop_body( idx[k] );
   }

   RAJA_FT_END ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over indirection array,
 *         including index count.
 *
 *         NOTE: lambda loop body requires two args (icount, index).
 *
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall_Icount(cilk_for_exec,
                   const Index_type* __restrict__ idx, const Index_type len,
                   const Index_type icount,
                   LOOP_BODY loop_body)
{
   RAJA_FT_BEGIN ;

   cilk_for ( Index_type k = 0 ; k < len ; ++k ) {
      loop_body( k+icount, idx[k] );
   }

   RAJA_FT_END ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over unstructured index set object.
 *
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall(cilk_for_exec,
            const ListSegment& is,
            LOOP_BODY loop_body)
{
   const Index_type* __restrict__ idx = is.getIndex();
   const Index_type len = is.getLength();

   RAJA_FT_BEGIN ;

   cilk_for ( Index_type k = 0 ; k < len ; ++k ) {
      loop_body( idx[k] );
   }

   RAJA_FT_END ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over unstructured index set object,
 *         including index count.
 *
 *         NOTE: lambda loop body requires two args (icount, index).
 *
 ******************************************************************************
 */
template <typename LOOP_BODY>
RAJA_INLINE
void forall_Icount(cilk_for_exec,
                   const ListSegment& is,
                   const Index_type icount,
                   LOOP_BODY loop_body)
{
   const Index_type* __restrict__ idx = is.getIndex();
   const Index_type len = is.getLength();

   RAJA_FT_BEGIN ;

   cilk_for ( Index_type k = 0 ; k < len ; ++k ) {
      loop_body( k+icount, idx[k] );
   }

   RAJA_FT_END ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  minloc reduction over given indirection array.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_minloc(cilk_for_exec,
                   const Index_type* __restrict__ idx, const Index_type len,
                   T* min, Index_type *loc,
                   LOOP_BODY loop_body)
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> min_tmp(nworkers);
   RAJAVec<Index_type> loc_tmp(nworkers);

   RAJA_FT_BEGIN ;

   for ( int i = 0; i < nworkers; ++i ) {
       min_tmp[i] = *min ;
       loc_tmp[i] = *loc ;
   }

   cilk_for ( Index_type k = 0 ; k < len ; ++k ) {
      loop_body( idx[k], &min_tmp[__cilkrts_get_worker_number()],
                         &loc_tmp[__cilkrts_get_worker_number()] );
   }

   for ( int i = 1; i < nworkers; ++i ) {
      if ( min_tmp[i] < min_tmp[0] ) {
         min_tmp[0] = min_tmp[i];
         loc_tmp[0] = loc_tmp[i];
      }
   }

   RAJA_FT_END ;

   *min = min_tmp[0] ;
   *loc = loc_tmp[0] ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  minloc reduction over unstructured index set object.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_minloc(cilk_for_exec,
                   const ListSegment& is,
                   T* min, Index_type *loc,
                   LOOP_BODY loop_body)
{
   forall_minloc(cilk_for_exec(),
                 is.getIndex(), is.getLength(),
                 min, loc,
                 loop_body);
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  maxloc reduction over given indirection array.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_maxloc(cilk_for_exec,
                   const Index_type* __restrict__ idx, const Index_type len,
                   T* max, Index_type *loc,
                   LOOP_BODY loop_body)
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> max_tmp(nworkers);
   RAJAVec<Index_type> loc_tmp(nworkers);

   RAJA_FT_BEGIN ;

   for ( int i = 0; i < nworkers; ++i ) {
       max_tmp[i] = *max ;
       loc_tmp[i] = *loc ;
   }

   cilk_for ( Index_type k = 0 ; k < len ; ++k ) {
      loop_body( idx[k], &max_tmp[__cilkrts_get_worker_number()],
                         &loc_tmp[__cilkrts_get_worker_number()] );
   }

   for ( int i = 1; i < nworkers; ++i ) {
      if ( max_tmp[i] > max_tmp[0] ) {
         max_tmp[0] = max_tmp[i];
         loc_tmp[0] = loc_tmp[i];
      }
   }

   RAJA_FT_END ;

   *max = max_tmp[0] ;
   *loc = loc_tmp[0] ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  maxloc reduction over unstructured index set object.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_maxloc(cilk_for_exec,
                   const ListSegment& is,
                   T* max, Index_type *loc,
                   LOOP_BODY loop_body)
{
   forall_maxloc(cilk_for_exec(),
                 is.getIndex(), is.getLength(),
                 max, loc,
                 loop_body);
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  sum reduction over given indirection array.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_sum(cilk_for_exec,
                const Index_type* __restrict__ idx, const Index_type len,
                T* sum,
                LOOP_BODY loop_body)
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> sum_tmp(nworkers);

   RAJA_FT_BEGIN ;

   for ( int i = 0; i < nworkers; ++i ) {
      sum_tmp[i] = 0 ;
   }

   cilk_for ( Index_type k = 0 ; k < len ; ++k ) {
      loop_body( idx[k], &sum_tmp[__cilkrts_get_worker_number()] );
   }

   RAJA_FT_END ;

   for ( int i = 0; i < nworkers; ++i ) {
      *sum += sum_tmp[i];
   }
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for  sum reduction over unstructured index set object.
 *
 ******************************************************************************
 */
template <typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_sum(cilk_for_exec,
                const ListSegment& is,
                T* sum,
                LOOP_BODY loop_body)
{
   forall_sum(cilk_for_exec(),
              is.getIndex(), is.getLength(),
              sum,
              loop_body);
}


//
//////////////////////////////////////////////////////////////////////
//
// The following function templates iterate over hybrid index set
// segments using cilk_for. Segment execution is defined by segment
// execution policy template parameter.
//
//////////////////////////////////////////////////////////////////////
//

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over segments of hybrid index set and 
 *         use execution policy template parameter to execute segments.
 *
 ******************************************************************************
 */
template <typename SEG_EXEC_POLICY_T,
          typename LOOP_BODY>
RAJA_INLINE
void forall( IndexSet::ExecPolicy<cilk_for_segit, SEG_EXEC_POLICY_T>,
             const IndexSet& is, LOOP_BODY loop_body )
{
   const int num_seg = is.getNumSegments();
   cilk_for ( int isi = 0; isi < num_seg; ++isi ) {
      SegmentType segtype = is.getSegmentType(isi);
      const void* iset = is.getSegmentISet(isi);

      switch ( segtype ) {

         case _RangeSeg_ : {
            forall(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const RangeSegment*>(iset)),
               loop_body
            );
            break;
         }

#if 0  // RDH RETHINK
         case _RangeStrideSeg_ : {
            forall(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const RangeStrideSegment*>(iset)),
               loop_body
            );
            break;
         }
#endif

         case _ListSeg_ : {
            forall(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const ListSegment*>(iset)),
               loop_body
            );
            break;
         }

         default : {
         }

      }  // switch on segment type

   } // iterate over parts of hybrid index set
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for iteration over segments of hybrid index set and
 *         use execution policy template parameter to execute segments.
 *
 *         This method passes count segment iteration.
 *
 *         NOTE: lambda loop body requires two args (icount, index).
 *
 ******************************************************************************
 */
template <typename SEG_EXEC_POLICY_T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_Icount( IndexSet::ExecPolicy<cilk_for_segit, SEG_EXEC_POLICY_T>,
                    const IndexSet& is, LOOP_BODY loop_body )
{
   const int num_seg = is.getNumSegments();
   cilk_for ( int isi = 0; isi < num_seg; ++isi ) {
      SegmentType segtype = is.getSegmentType(isi);
      const void* iset = is.getSegmentISet(isi);
      Index_type icount = is.getSegmentIcount(isi);

      switch ( segtype ) {

         case _RangeSeg_ : {
            foral_Icount(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const RangeSegment*>(iset)),
               icount,
               loop_body
            );
            break;
         }

#if 0  // RDH RETHINK
         case _RangeStrideSeg_ : {
            foral_Icount(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const RangeStrideSegment*>(iset)),
               icount,
               loop_body
            );
            break;
         }
#endif

         case _ListSeg_ : {
            foral_Icount(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const ListSegment*>(iset)),
               icount,
               loop_body
            );
            break;
         }

         default : {
         }

      }  // switch on segment type

   } // iterate over parts of hybrid index set
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for minloc reduction over segments of hybrid index set and
 *         use execution policy template parameter to execute segments.
 *
 ******************************************************************************
 */
template <typename SEG_EXEC_POLICY_T,
          typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_minloc( IndexSet::ExecPolicy<cilk_for_segit, SEG_EXEC_POLICY_T>,
                    const IndexSet& is,
                    T* min, Index_type *loc,
                    LOOP_BODY loop_body )
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> min_tmp(nworkers);
   RAJAVec<Index_type> loc_tmp(nworkers);

   for ( int i = 0; i < nworkers; ++i ) {
       min_tmp[i] = *min ;
       loc_tmp[i] = *loc ;
   }

   const int num_seg = is.getNumSegments();
   cilk_for ( int isi = 0; isi < num_seg; ++isi ) {
      SegmentType segtype = is.getSegmentType(isi);
      const void* iset = is.getSegmentISet(isi);

      switch ( segtype ) {

         case _RangeSeg_ : {
            forall_minloc(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const RangeSegment*>(iset)),
               &min_tmp[__cilkrts_get_worker_number()],
               &loc_tmp[__cilkrts_get_worker_number()],
               loop_body
            );
            break;
         }

#if 0  // RDH RETHINK
         case _RangeStrideSeg_ : {
            forall_minloc(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const RangeStrideSegment*>(iset)),
               &min_tmp[__cilkrts_get_worker_number()],
               &loc_tmp[__cilkrts_get_worker_number()],
               loop_body
            );
            break;
         }
#endif

         case _ListSeg_ : {
            forall_minloc(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const ListSegment*>(iset)),
               &min_tmp[__cilkrts_get_worker_number()],
               &loc_tmp[__cilkrts_get_worker_number()],
               loop_body
            );
            break;
         }

         default : {
         }

      }  // switch on segment type

   } // iterate over segments of hybrid index set


   for ( int i = 1; i < nworkers; ++i ) {
      if ( min_tmp[i] < min_tmp[0] ) {
         min_tmp[0] = min_tmp[i];
         loc_tmp[0] = loc_tmp[i];
      }
   }

   *min = min_tmp[0] ;
   *loc = loc_tmp[0] ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for maxloc  reduction over segments of hybrid index set and
 *         use execution policy template parameter to execute segments.
 *
 ******************************************************************************
 */
template <typename SEG_EXEC_POLICY_T,
          typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_maxloc( IndexSet::ExecPolicy<cilk_for_segit, SEG_EXEC_POLICY_T>,
                    const IndexSet& is,
                    T* max, Index_type *loc,
                    LOOP_BODY loop_body )
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> max_tmp(nworkers);
   RAJAVec<Index_type> loc_tmp(nworkers);

   for ( int i = 0; i < nworkers; ++i ) {
       max_tmp[i] = *max ;
       loc_tmp[i] = *loc ;
   }

   const int num_seg = is.getNumSegments();
   cilk_for ( int isi = 0; isi < num_seg; ++isi ) {
      SegmentType segtype = is.getSegmentType(isi);
      const void* iset = is.getSegmentISet(isi);

      switch ( segtype ) {

         case _RangeSeg_ : {
            forall_maxloc(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const RangeSegment*>(iset)),
               &max_tmp[__cilkrts_get_worker_number()],
               &loc_tmp[__cilkrts_get_worker_number()],
               loop_body
            );
            break;
         }

#if 0  // RDH RETHINK
         case _RangeStrideSeg_ : {
            forall_maxloc(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const RangeStrideSegment*>(iset)),
               &max_tmp[__cilkrts_get_worker_number()],
               &loc_tmp[__cilkrts_get_worker_number()],
               loop_body
            );
            break;
         }
#endif

         case _ListSeg_ : {
            forall_maxloc(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const ListSegment*>(iset)),
               &max_tmp[__cilkrts_get_worker_number()],
               &loc_tmp[__cilkrts_get_worker_number()],
               loop_body
            );
            break;
         }

         default : {
         }

      }  // switch on segment type

   } // iterate over segments of hybrid index set

   for ( int i = 1; i < nworkers; ++i ) {
      if ( max_tmp[i] > max_tmp[0] ) {
         max_tmp[0] = max_tmp[i];
         loc_tmp[0] = loc_tmp[i];
      }
   }

   *max = max_tmp[0] ;
   *loc = loc_tmp[0] ;
}

/*!
 ******************************************************************************
 *
 * \brief  cilk_for sum  reduction over segments of hybrid index set and
 *         use execution policy template parameter to execute segments.
 *
 ******************************************************************************
 */
template <typename SEG_EXEC_POLICY_T,
          typename T,
          typename LOOP_BODY>
RAJA_INLINE
void forall_sum( IndexSet::ExecPolicy<cilk_for_segit, SEG_EXEC_POLICY_T>,
                 const IndexSet& is,
                 T* sum,
                 LOOP_BODY loop_body )
{
   const int nworkers = __cilkrts_get_nworkers();

   RAJAVec<T> sum_tmp(nworkers);

   for ( int i = 0; i < nworkers; ++i ) {
      sum_tmp[i] = 0 ;
   }

   const int num_seg = is.getNumSegments();
   cilk_for ( int isi = 0; isi < num_seg; ++isi ) {
      SegmentType segtype = is.getSegmentType(isi);
      const void* iset = is.getSegmentISet(isi);

      switch ( segtype ) {

         case _RangeSeg_ : {
            forall_sum(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const RangeSegment*>(iset)),
               &sum_tmp[__cilkrts_get_worker_number()],
               loop_body
            );
            break;
         }

#if 0  // RDH RETHINK
         case _RangeStrideSeg_ : {
            forall_sum(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const RangeStrideSegment*>(iset)),
               &sum_tmp[__cilkrts_get_worker_number()],
               loop_body
            );
            break;
         }
#endif

         case _ListSeg_ : {
            forall_sum(
               SEG_EXEC_POLICY_T(),
               *(static_cast<const ListSegment*>(iset)),
               &sum_tmp[__cilkrts_get_worker_number()],
               loop_body
            );
            break;
         }

         default : {
         }

      }  // switch on segment type

   } // iterate over segments of hybrid index set

   for ( int i = 0; i < nworkers; ++i ) {
      *sum += sum_tmp[i];
   }
}


}  // closing brace for RAJA namespace

#endif  // closing endif for header file include guard

