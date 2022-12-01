#include "RAJA_test-base.hpp"

template <typename EXEC_POLICY>
void SymSegmentTestImpl() {

  EXPECT_EQ(1,0);
}

TYPED_TEST_SUITE_P(SymSegmentTest);
template <typename T>
class SymSegmentTest : public ::testing::Test
{
};

TYPED_TEST_P(SymSegmentTest, SymSegment)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  SymSegmentTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(SymSegmentTest, SymSegment);
