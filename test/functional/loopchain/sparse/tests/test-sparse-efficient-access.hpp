#include "RAJA_test-base.hpp"
#include <string>
#include <vector>


template <typename EXEC_POLICY>
void SparseEfficientAccessTestImpl() {
  
  using namespace RAJA;
  idx_t N = 10;
  auto A = make_random_sparse_view2<double>(N);
  double s01 = 0.0;
  double s10 = 0.0;
  using POL = KernelPolicy<statement::For<0,loop_exec, statement::For<1,loop_exec, statement::Lambda<0>>>>;

  auto dense_space = make_tuple(RangeSegment(0,N), RangeSegment(0,N));
  std::vector<idx_t> perm01{0,1};
  std::vector<idx_t> perm10{1,0};
  auto sparse_space01 = sparsify_iteration_space<POL>(dense_space, A, perm01);
  auto sparse_space10 = sparsify_iteration_space<POL>(dense_space, A, perm10);
  A.set_permutation(perm01);
  
  auto lam01 = [&](auto i, auto j) {
    s01 += A(i,j);
  };
  auto lam10 = [&](auto i, auto j) {
    s10 += A(j,i);
  };

  kernel<POL>(sparse_space01, lam01);

  EXPECT_EQ(A.get_hit_rate(), 1.0);
  A.reset_counters();
  
  kernel<POL>(sparse_space10, lam10);
  EXPECT_NE(A.get_hit_rate(), 1.0);
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
}

REGISTER_TYPED_TEST_SUITE_P(SparseEfficientAccessTest, SparseEfficientAccess);
