#include "RAJA_test-base.hpp"


template <typename EXEC_POLICY>
void SymBoundTestImpl() {

  using namespace RAJA;
  idx_t Ni = 10;
  idx_t Nj = 10;
  
  
  auto lo = make_symbolic_bound(0);
  auto hi = make_symbolic_bound(Ni);
  
  EXPECT_EQ(lo->value(), 0);
  EXPECT_EQ(hi->value(), Ni);
  
  auto i_seg = make_symbolic_segment(0,Ni);
  auto segbound = make_symbolic_bound(i_seg); 

  EXPECT_EQ(segbound->value(), 0);
  i_seg.currVal ++;
  EXPECT_EQ(i_seg.currVal, 1);
  EXPECT_EQ(segbound->value(), i_seg.currVal);
  i_seg.currVal --;
 
  auto i_seg2 = i_seg * 2;
  auto seg2bound = make_symbolic_bound(i_seg2);
  EXPECT_EQ(seg2bound->value(), 0);
  i_seg.currVal ++;
  EXPECT_EQ(seg2bound->value(), 2);

}

template <typename EXEC_POLICY>
void SymSegmentTestImpl() {

  using namespace RAJA;
  idx_t Ni = 10;
  idx_t Nj = 10;
  
  auto i_seg = make_symbolic_segment(0,Ni);
  auto j1_seg = make_symbolic_segment(i_seg,Nj);
  auto j2_seg = make_symbolic_segment(i_seg, i_seg * 2);
  
  std::stringstream s;
  auto lambda1 = [&](auto i) {
    s << i;
  };

  //auto forall1 = make_forall<seq_exec>(i_seg, lambda1);
  
  //forall1();
  //EXPECT_EQ(s.str().c_str(), "0123456789");



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
  SymBoundTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(SymSegmentTest, SymSegment);
