#ifndef _RAJA_FDBLAS_HPP
#define _RAJA_FDBLAS_HPP

#include "RAJA/config.hpp"

namespace RAJA {


// creates a kernel that computes C = alpha * A * B + beta * C
template <typename T>
auto gemm(View<T,Layout<2>> &A, View<T,Layout<2>> &B, View<T,Layout<2>> &C, T alpha, T beta) {
  auto Ni = C.get_layout().sizes[0];
  auto Nj = C.get_layout().sizes[1];
  auto Nk = A.get_layout().sizes[0];
  
  auto segments = make_tuple(RangeSegment(0,Ni), RangeSegment(0,Nj), RangeSegment(0,Nk));
  
  using POLICY = KernelPolicy<
    statement::For<0,omp_parallel_for_exec,
      statement::For<1,seq_exec,
        statement::Lambda<0>,
        statement::For<2,simd_exec,
          statement::Lambda<1>
  >>>>;
  auto initC = [&](auto i, auto j) {
    C(i,j) *= beta;
  };
  auto comp = [&](auto i, auto j, auto k) {
    C(i,j) += A(i,k) * B(k,j) * alpha;
  };

  return make_kernel<POLICY>(segments, initC, comp);
}

//computes y = alpha * A * x + beta * y
template <typename T>
auto gemv(View<T,Layout<2>> &A, View<T,Layout<1>> &x, View<T,Layout<1>> y, T alpha, T beta) {

  auto yInit = [&](auto i) {
    y(i) *= beta;
  };
  auto comp = [&](auto i, auto j) {
    y(i) += A(i,j) * x(j) * alpha;
  };

  using POLICY = KernelPolicy<
    statement::For<0,omp_parallel_for_exec,
      statement::Lambda<0>,
      statement::For<1,simd_exec,
        statement::Lambda<1>
      >
    >
  >;

  auto Ni = A.get_layout().sizes[0];
  auto Nj = A.get_layout().sizes[1];
  auto segments = make_tuple(RangeSegment(0,Ni), RangeSegment(0,Nj));
  return make_kernel<POLICY>(segments, yInit, comp);
}

//computes A = alpha * x * y^T + A
template <typename T>
auto ger(View<T,Layout<2>> &A, View<T,Layout<1>> &x, View<T,Layout<1>> &y, T alpha) {
  auto Ni = A.get_layout().sizes[0];
  auto Nj = A.get_layout().sizes[1];
  auto segments = make_tuple(RangeSegment(0,Ni), RangeSegment(0,Nj));
  
  using POLICY = KernelPolicy<
    statement::For<0,omp_parallel_for_exec,
      statement::For<1,simd_exec,
        statement::Lambda<0>
  >>>;

  auto comp = [&](auto i, auto j) {
    A(i,j) += alpha * x(i) * y(j);
  };

  return make_kernel<POLICY>(segments, comp);
}

}

#endif
