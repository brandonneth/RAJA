#include "RAJA_test-base.hpp"

template <typename EXEC_POLICY>
void MortonTestImpl() {

  using namespace RAJA;
  size_t dimN = 4;
  size_t tileN = 1;
  std::array<size_t, 2> dimSizes{{dimN,dimN}}; 
  std::array<size_t, 2> tileSizes{{tileN,tileN}};
  
  auto mortonLayout = make_morton_layout<2>(dimSizes, tileSizes);

  EXPECT_EQ(mortonLayout(0,0), 0);
  EXPECT_EQ(mortonLayout(0,1), 1);
  EXPECT_EQ(mortonLayout(1,0), 2);
  EXPECT_EQ(mortonLayout(1,1), 3);
  EXPECT_EQ(mortonLayout(0,2), 4);
  EXPECT_EQ(mortonLayout(0,3), 5);
  EXPECT_EQ(mortonLayout(1,2), 6);
  EXPECT_EQ(mortonLayout(1,3), 7);
  EXPECT_EQ(mortonLayout(2,0), 8);
  EXPECT_EQ(mortonLayout(2,1), 9);
  EXPECT_EQ(mortonLayout(3,0), 10);
  EXPECT_EQ(mortonLayout(3,1), 11);
  EXPECT_EQ(mortonLayout(2,2), 12);
  EXPECT_EQ(mortonLayout(2,3), 13);
  EXPECT_EQ(mortonLayout(3,2), 14);
  EXPECT_EQ(mortonLayout(3,3), 15);

  auto tiledMortonLayout = make_morton_layout<2,size_t>({{8,8}}, {{4,4}});

  EXPECT_EQ(tiledMortonLayout(0,0), 0); 
  EXPECT_EQ(tiledMortonLayout(0,3), 3); 
  EXPECT_EQ(tiledMortonLayout(1,3), 7); 
  EXPECT_EQ(tiledMortonLayout(3,3), 15); 
  EXPECT_EQ(tiledMortonLayout(0,4), 16); 
  EXPECT_EQ(tiledMortonLayout(0,6), 18); 
  EXPECT_EQ(tiledMortonLayout(2,6), 26); 
  EXPECT_EQ(tiledMortonLayout(4,0), 32); 
  EXPECT_EQ(tiledMortonLayout(4,1), 33); 
  EXPECT_EQ(tiledMortonLayout(4,4), 48); 
  EXPECT_EQ(tiledMortonLayout(7,7), 63); 

}

TYPED_TEST_SUITE_P(MortonTest);
template <typename T>
class MortonTest : public ::testing::Test
{
};

TYPED_TEST_P(MortonTest, Morton)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  MortonTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(MortonTest, Morton);
