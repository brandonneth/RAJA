#include "RAJA_test-base.hpp"

template <typename EXEC_POLICY>
void FormatDecisionsTestImpl() {
  
  using namespace RAJA;
  double * _a = new double[9];
  for(int i = 0; i < 9; i++) {
    _a[i] = i;
  }

  using VIEW2 = RAJA::View<double, RAJA::Layout<2,RAJA::idx_t>>;

 
  

  std::array<RAJA::idx_t, 2> perm_01 {{0,1}};
  std::array<RAJA::idx_t, 2> perm_10 {{1,0}};
  auto layout_01 = RAJA::make_permuted_layout( {{3,3}}, perm_01 );
  auto layout_10 = RAJA::make_permuted_layout( {{3,3}}, perm_10 );
  VIEW2 a(_a,layout_01);
  permute_view(a, layout_10)();

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      EXPECT_EQ(a(i,j), i*3 + j);
      auto val = j*3 + i;
      EXPECT_EQ(_a[i*3+j], val);
    }
  }

  

}

TYPED_TEST_SUITE_P(FormatDecisionsTest);
template <typename T>
class FormatDecisionsTest : public ::testing::Test
{
};

TYPED_TEST_P(FormatDecisionsTest, FormatDecisions)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  FormatDecisionsTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(FormatDecisionsTest, FormatDecisions);
