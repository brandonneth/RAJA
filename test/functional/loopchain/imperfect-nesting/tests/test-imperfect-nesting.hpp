#include "RAJA_test-base.hpp"

#include <unordered_set>
template <typename EXEC_POLICY>
void ImperfectNestingTestImpl() {
  
  using namespace RAJA;
  
  size_t N = 8;


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
  auto lam1 = [&](auto i, auto) {
    s << i;
  };
  auto lam2 = [&](auto, auto j) {
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
      statement::Lambda<0, RAJA::Segs<0>>,
      statement::For<1,loop_exec,
        statement::Lambda<1, RAJA::Segs<0,1>>
      >
    >
  >;

  for(size_t i = 0; i < N; i++) {
    correct << i;
    for(size_t j = 0; j < N; j++) {
      correct << j;
    }
  }

  auto lam1_1 = [&](auto i) {
    s << i;
  };
  auto knl2 = make_kernel<forlamforlam>(segs, lam1_1, lam2);
  knl2();

  EXPECT_STREQ(correct.str().c_str(), s.str().c_str());



  
  using VIEW2 = RAJA::View<double, RAJA::Layout<2,RAJA::idx_t>>;
  VIEW2 a(new double[N*N], N, N);
  VIEW2 b(new double[N*N], N, N);
  VIEW2 c(new double[N*N], N, N);

  auto symlam1 = [&](auto i, auto j) {
    a(i,j) = 0;
    b(i,j) = 1;
  };
  auto symlam2 = [&](auto i, auto j) {
    c(i,j) = 2;
  };

  auto symknl1 = RAJA::make_kernel<forforlamlam>(segs, symlam1, symlam2);
  
  auto accesses1 = symknl1.execute_symbolically();
  EXPECT_EQ(accesses1.size(), 3) << "forforlamlam symbolic evaluation does not evaluate both lambdas.";

  for(auto access: accesses1) {
    if (access.is_view(a)) {
      EXPECT_EQ(access.numIterations, N*N) << "symknl1's access to a counts iterations wrong.";
    } else if (access.is_view(b)) {
      EXPECT_EQ(access.numIterations, N*N) << "symknl1's access to b counts iterations wrong.";
    } else if (access.is_view(c)) {
      EXPECT_EQ(access.numIterations, N*N) << "symknl1's access to c counts iterations wrong";
    } else {
      EXPECT_EQ(0,1) << "Access object's view pointer does not point to any of the views.";
    } 
  }

  auto symlam1_1 = [&](auto i) {
    a(i,0) = 0;
  };


  
  auto symknl2 = RAJA::make_kernel<forlamforlam>(segs, symlam1_1, symlam2);

  auto accesses2 = symknl2.execute_symbolically();
  EXPECT_EQ(accesses2.size(), 2) << "forlamforlam symbolic evaluation does not evaluate both lambdas.";
  for(auto access: accesses2) {
    if (access.is_view(a)) {
      EXPECT_EQ(access.numIterations, N) << "symknl2's access to a counts iterations wrong.";
    } else if (access.is_view(c)) {
      EXPECT_EQ(access.numIterations, N*N) << "symknl2's access to c counts iterations wrong";
    } else {
      EXPECT_EQ(0,1) << "Access object's view pointer does not point to any of the views.";
    } 
  }

  

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
