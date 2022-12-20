#include "RAJA_test-base.hpp"

#include <unordered_set>
template <typename EXEC_POLICY>
void ImperfectNestingTestImpl() {
  
  using namespace RAJA;
  double * _a = new double[9];

  for(int i = 0; i < 9; i++) {
    _a[i] = i;
  }

  size_t N = 8;
  using VIEW2 = RAJA::View<double, RAJA::Layout<2,RAJA::idx_t>>;


  auto segs = make_tuple(RangeSegment(0,N), RangeSegment(0,N));

  using forforlamlam = KernelPolicy<
    statement::For<0,loop_exec,
      statement::For<1,loop_exec,
        statement::Lambda<0>,
        statement::Lambda<1>
      >
    >
  >;
  std::stringstream s;
  auto lam1 = [&](auto i, auto j) {
    s << i;
  };
  auto lam2 = [&](auto i, auto j) {
    s << j;
  };
  std::stringstream correct;
  for(size_t i = 0; i < N; i++) {
    for(size_t j = 0; j  <N; j++) {
      correct << i;
      correct << j;
    }
  }

  auto knl1 = make_kernel<forforlamlam>(segs, lam1, lam2);

  knl1();

  EXPECT_STREQ(correct.str().c_str(), s.str().c_str());
 
  correct.str("");
  correct.clear();
  s.str("");
  s.clear();
  using forlamforlam = KernelPolicy<
    statement::For<0,loop_exec,
      statement::Lambda<0>,
      statement::For<1,loop_exec,
        statement::Lambda<1>
      >
    >
  >;

  for(size_t i = 0; i < N; i++) {
    correct << i;
    for(size_t j = 0; j < N; j++) {
      correct << j;
    }
  }

  auto knl2 = make_kernel<forlamforlam>(segs, lam1, lam2);
  knl2();
  EXPECT_STREQ(correct.str().c_str(), s.str().c_str());

}

TYPED_TEST_SUITE_P(ImperfectNestingTest);
template <typename T>
class ImperfectNestingTest : public ::testing::Test
{
};

TYPED_TEST_P(ImperfectNestingTest, ImperfectNesting)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  ImperfectNestingTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(ImperfectNestingTest, ImperfectNesting);
