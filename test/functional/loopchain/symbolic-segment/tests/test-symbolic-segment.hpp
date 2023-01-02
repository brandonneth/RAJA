#include "RAJA_test-base.hpp"


template <typename EXEC_POLICY>
void SymBoundTestImpl() {

  using namespace RAJA;
  idx_t Ni = 10;
  
  
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


  using VIEW2 = RAJA::View<double, RAJA::Layout<2,RAJA::idx_t>>;
  idx_t N = 10;
  VIEW2 a(new double[N*N], N,N);
  VIEW2 b(new double[N*N], N,N);


  auto seg1 = make_symbolic_segment(0,N,"seg1");
  auto seg2 = make_symbolic_segment(seg1, N, "seg2");
  
  auto parlam = [&](auto i, auto j) {
    b(i,j) = a(i,j);
  };

  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      a(i,j) = 1;
      b(i,j) = 2;
    }
  }

  using pol2d_par = KernelPolicy<
   statement::For<0,omp_parallel_for_exec,
      statement::For<1,loop_exec,
        statement::Lambda<0>
      >
    >
  >;

  auto par_knl = make_kernel<pol2d_par>(make_tuple(seg1,seg2), parlam);

  par_knl();

  for(int i = 0; i < N; i++) {
    for(int j = 0; j < i; j++) {
      EXPECT_EQ(b(i,j), 2) << "B view assigned on iteration " << i << ", " << j ;
    }
    for(int j = i; j < N; j++) {
      EXPECT_EQ(b(i,j), 1) << "B view not assigned on iteration " << i << ", " << j ;
    }
  }

  auto accesses = par_knl.execute_symbolically();

  size_t numIterations =N*(N+1) / 2;
  for (auto access : accesses) {
    if (access.is_view(a)) {
      EXPECT_EQ(access.numIterations, numIterations) << "Access to a counted iterations wrong.";
    }
    else if (access.is_view(b)) {
      EXPECT_EQ(access.numIterations, numIterations) << "Access to b counted iterations wrong.";
    } else {
      EXPECT_EQ(0,1) << "Access object's view reference is not a or b.";
    }
  }
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
