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
  i_seg.increment();
  EXPECT_EQ(i_seg.value(), 1);
  EXPECT_EQ(segbound->value(), i_seg.value());
  i_seg.decrement();
 
  auto i_seg2 = i_seg * 2;
  auto seg2bound = make_symbolic_bound(i_seg2);
  EXPECT_EQ(seg2bound->value(), 0);
  i_seg.increment();
  EXPECT_EQ(seg2bound->value(), 2);

}

template <typename EXEC_POLICY>
void SymSegmentTestImpl() {

  using namespace RAJA;
  idx_t Ni = 3;
  idx_t Nj = 3;
  
  auto i_seg = make_symbolic_segment(0,Ni,"i_seg");
  auto j1_seg = make_symbolic_segment(i_seg,Nj, "j1_seg");
  auto j2_seg = make_symbolic_segment(i_seg, i_seg * 2, "j2_seg");
  
  std::stringstream s;
  auto lambda1 = [&](auto i) {
    s << i;
  };

  auto forall1 = make_forall<seq_exec>(i_seg, lambda1);
  
  forall1();
  EXPECT_STREQ(s.str().c_str(), "012");
  i_seg.reset();
  using pol2d = KernelPolicy<
    statement::For<0,loop_exec,
      statement::For<1,loop_exec,
        statement::Lambda<0>
      >
    >
  >;

  s.str(std::string());
  auto knl1 = make_kernel<pol2d>(make_tuple(i_seg, RangeSegment(0,3)), [&](auto i, auto j) {
    s << i << j;
  });
  knl1();
  EXPECT_STREQ(s.str().c_str(), "000102101112202122");
  i_seg.reset();
  s.str(std::string());

  EXPECT_EQ(i_seg.value(), 0) << "i segment did not reset\n";
  EXPECT_EQ(j1_seg.value(), 0) << "j segment does not have correct starting value before its first use.\n";




  auto knl2 = make_kernel<pol2d>(make_tuple(i_seg, j1_seg), [&](auto i, auto j) {
    s << i << j;
  });
  std::cout << "Launching kernel with two symbolic segments\n";
  knl2();
  EXPECT_STREQ(s.str().c_str(), "000102111222");
  i_seg.reset();
  j1_seg.reset();
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
