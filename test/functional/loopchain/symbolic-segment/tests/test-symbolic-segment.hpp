#include "RAJA_test-base.hpp"

template <typename EXEC_POLICY>
void SymbolicSegmentTestImpl() {

  EXPECT_EQ(1,1);
}

TYPED_TEST_SUITE_P(SymbolicSegmentTest);
template <typename T>
class SymbolicSegmentTest : public ::testing::Test
{
};

TYPED_TEST_P(SymbolicSegmentTest, SymbolicSegment)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  SymbolicSegmentTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(SymbolicSegmentTest, SymbolicSegment);
