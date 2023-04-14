#include "RAJA_test-base.hpp"
#include <string>
#include <vector>


template <typename EXEC_POLICY>
void TestTrimSegmentTuple() {
  using namespace RAJA;

  auto seg1 = RangeSegment(0,10);
  auto seg2 = RangeSegment(10,20);
  auto seg3 = RangeSegment(20,30);


  std::vector<idx_t> p10{1,0};
  std::vector<idx_t> p01{0,1};
  std::vector<idx_t> p20{2,0};
  std::vector<idx_t> p12{1,2};


  auto segs1 = camp::tie(seg1,seg2,seg3);
  auto trimmed1_01 = trim_segment_tuple<2>(segs1, p01);
  EXPECT_EQ(camp::get<0>(trimmed1_01).begin()[0], 0);
  EXPECT_EQ(camp::get<0>(trimmed1_01).end()[0], 10);
  EXPECT_EQ(camp::get<1>(trimmed1_01).begin()[0], 10);
  EXPECT_EQ(camp::get<1>(trimmed1_01).end()[0], 20);

  auto trimmed1_10 = trim_segment_tuple<2>(segs1, p10);
  
  EXPECT_EQ(camp::get<0>(trimmed1_01).begin()[0], 0);
  EXPECT_EQ(camp::get<0>(trimmed1_01).end()[0], 10);
  EXPECT_EQ(camp::get<1>(trimmed1_01).begin()[0], 10);
  EXPECT_EQ(camp::get<1>(trimmed1_01).end()[0], 20);

  
  

}

TYPED_TEST_SUITE_P(SparseSparsifyUnitTest);
template <typename T>
class SparseSparsifyUnitTest : public ::testing::Test
{
};

TYPED_TEST_P(SparseSparsifyUnitTest, SparseSparsifyUnit)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  TestTrimSegmentTuple<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SparseSparsifyUnitTest, SparseSparsifyUnit);
