#include "RAJA_test-base.hpp"
#include <string>
#include <vector>


template <typename EXEC_POLICY>
void GauSeiDIAGMisses() {
  using namespace RAJA;
  using DenseView1 = View<double, Layout<1>>;
  idx_t dimSize = 10;
  double nonzeroDensity = 1.0;

  auto refData = make_random_sparse_view<double, 2, DIAG>(dimSize, nonzeroDensity, true);
  DenseView1 x(new double[dimSize], dimSize);
  DenseView1 temp(new double[dimSize], dimSize);
  DenseView1 b(new double[dimSize], dimSize);

  std::cout << "running gausie diagonal misses.";

  for(int i = 0; i < dimSize; i++) {
    //EXPECT_NE(refData(i,i), 0) << "zero on diagonal.\n";
    x(i) = 1.0;
    temp(i) = 0;
    b(i) = 1.0;
  }
  using POLICY = KernelPolicy<
      statement::For<0,loop_exec,
        statement::Lambda<0,Segs<0>>,
        statement::For<1,loop_exec,
          statement::Lambda<1>
        >,
        statement::Lambda<2,Segs<0>>
      >
    >;
  
  auto seg1 = RangeSegment(0,dimSize);
  auto seg2 = RangeSegment(0,dimSize);
  auto dense_segs = make_tuple(seg1, seg2);

    auto lam1 = [&](auto i) {
      temp(i) = 0.0;
    };
    auto lam2 = [&](auto i, auto j) {
      if (j != i) {
        auto summand = refData(i,j) * x(j);
        double newVal = temp(i) + summand;
        temp(i) = newVal;
      }
    };
    auto lam3 = [&](auto i) {
      x(i) = (b(i) - temp(i)) / refData(i,i);
    };

  auto knl = make_sparse_kernel_sym<POLICY, 1>(make_tuple(0,1), dense_segs, refData, lam1, lam2, lam3);
  knl();
  EXPECT_EQ(refData.get_misses(), 0);
  knl();
  EXPECT_EQ(refData.get_misses(), 0) << "fails on second one.";
  knl();
  EXPECT_EQ(refData.get_misses(), 0) << "fails on third one.";

}

template <typename EXEC_POLICY>
void SparseEfficientAccessTestImpl() {
  
  using namespace RAJA;
  idx_t N = 10;
  auto A = make_random_sparse_view<double,2,COO>(N,0.1);
  double s01 = 0.0;
  double s10 = 0.0;
  using POL = KernelPolicy<statement::For<0,loop_exec, statement::For<1,loop_exec, statement::Lambda<0>>>>;

  auto dense_space = make_tuple(RangeSegment(0,N), RangeSegment(0,N));
  std::vector<idx_t> perm01{0,1};
  std::vector<idx_t> perm10{1,0};
  auto sparse_space01 = sparsify_iteration_space_compressed<POL>(dense_space, A, perm01);
  auto sparse_space10 = sparsify_iteration_space_compressed<POL>(dense_space, A, perm10);
  A.set_permutation(perm01);
  
  auto lam01 = [&](auto i, auto j) {
    s01 += A(i,j);
  };
  auto lam10 = [&](auto i, auto j) {
    s10 += A(j,i);
  };

  std::vector<idx_t> dim0 = A.dim(0);
  std::vector<idx_t> dim1 = A.dim(1);
  for(size_t i = 0; i < dim0.size(); i++) {
  }
  std::cout << "\n";
  for(size_t i = 0; i < dim0.size(); i++) {
  }

  kernel<POL>(sparse_space01, lam01);

  EXPECT_EQ(A.get_hit_rate(), 1.0) << "efficient access failed on access A(i,j)\n";
  A.reset_counters();
  
  kernel<POL>(sparse_space10, lam10);
  EXPECT_NE(A.get_hit_rate(), 1.0) << "efficient access failed on access A(j,i)\n";
  A.reset_counters(); 
}

TYPED_TEST_SUITE_P(SparseEfficientAccessTest);
template <typename T>
class SparseEfficientAccessTest : public ::testing::Test
{
};

TYPED_TEST_P(SparseEfficientAccessTest, SparseEfficientAccess)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  SparseEfficientAccessTestImpl<EXEC_POLICY>();
  GauSeiDIAGMisses<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SparseEfficientAccessTest, SparseEfficientAccess);
