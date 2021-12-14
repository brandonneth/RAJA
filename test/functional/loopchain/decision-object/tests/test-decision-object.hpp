#include "RAJA_test-base.hpp"

template <typename EXEC_POLICY>
void DecisionObjectTestImpl() {

  using namespace RAJA;
  RAJA::idx_t N = 16;
  double * _a = new double[N*N];
  double * _b = new double[N*N];
  double * _c = new double[N*N];
  double * _d = new double[N*N];
  double * _e = new double[N*N];

  using KPOL = KernelPolicy<
    statement::For<0, loop_exec,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
          statement::Lambda<0>
        >
      >
    >
  >;
  
  using VIEW2 = View<double, RAJA::Layout<2>>;
  VIEW2 a(_a, N, N);
  VIEW2 b(_b, N, N);
  VIEW2 c(_c, N, N);
  VIEW2 d(_d, N, N);
  VIEW2 e(_e, N, N);

  auto initialize = [&]() {
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      a(i,j) = 1;
      b(i,j) = 2;
      c(i,j) = 0;
      d(i,j) = 5;
      e(i,j) = 0;
    }
  }
  for(int i = 0; i < N/2; i++) {
    a(i,0) += 2;
    b(i,i) *= -1;
  }
  a(0,1) = -3;
  };
 
  auto checksum = [&]() {
    double s = 0;
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        s += i*c(i,j) + j*e(i,j);
      }
    }
    return s;
  }; 
  auto lambda1 = [&](auto i0, auto i1, auto i2) {
    c(i0,i2) += a(i0,i1) * b(i1, i2);
  };

  auto lambda2 = [&](auto i0, auto i1, auto i2) {
    e(i1,i0) += d(i0,i2) * c(i2,i0);
  };

   auto print_array = [&](auto a) {
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        std::cout << a(i,j) << " ";
      }
      std::cout << "\n";
    }

  };

  auto iterspace = make_tuple(RangeSegment(0,N), RangeSegment(0,N), RangeSegment(0,N));

  auto knl1 = make_kernel<KPOL>(iterspace, lambda1);

  auto knl2 = make_kernel<KPOL>(iterspace, lambda2);

  //tie() creates tuple of references
  //
  auto decisions = create_layout_decisions<2>(tie(c));

  decisions.add_decision(c, c.get_layout(), 0);
  
  std::array<camp::idx_t, 2> cs{c.get_layout().sizes[0],c.get_layout().sizes[1]};


  auto c_perm = make_permuted_layout(cs, {{1,0}});

  decisions.add_decision(c, c_perm, 1); 
  decisions.add_decision(c, c_perm, 2); 

  auto computation = manipulate_layouts(decisions, knl1, knl2);
  
   
  initialize();
  knl1();
  knl2();
  auto baseline = checksum();

  auto b_indexed = c.get_data()[1];
/* 
  std::cout << "After baseline computation.\n";
  std::cout << "View c:\n";
  print_array(c);
  std::cout << "View c's data memory, index 1: " << c.get_data()[1] << "\n";
*/
 
  initialize();
  computation();
  auto referent = checksum();

/*
  std::cout << "After manipulated computation.\n";
  std::cout << "View c:\n";
  print_array(c);
  std::cout << "View c's data memory, index 1: " << c.get_data()[1] << "\n";
*/
  auto r_indexed = c.get_data()[1];
  
  EXPECT_EQ(baseline, referent);
  EXPECT_NE(b_indexed, r_indexed) << "c array layout did not change between baseline and manipulated computations\n";
 
}

TYPED_TEST_SUITE_P(DecisionObjectTest);
template <typename T>
class DecisionObjectTest : public ::testing::Test
{
};

TYPED_TEST_P(DecisionObjectTest, DecisionObject)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  DecisionObjectTestImpl<EXEC_POLICY>();

}


TEST(DecisionObjectTest, InitializedWithLayouts) {
  using namespace RAJA;
  RAJA::idx_t N = 16;
  double * _a = new double[N*N];
 
  using KPOL = KernelPolicy<
    statement::For<0, loop_exec,
      statement::For<1, loop_exec,
        statement::Lambda<0>
      >
    >
  >;
  
  using VIEW2 = View<double, RAJA::Layout<2>>;
  
  std::array<camp::idx_t, 2> sizes{{N,N}};
  std::array<camp::idx_t, 2> perm{{1,0}};
  Layout<2> flipped = make_permuted_layout(sizes, perm);
  VIEW2 a(_a, flipped);

  auto l = [&](auto i, auto j) {
    a(i,j) = 1;
  };

  auto segment = make_tuple(RangeSegment(0,N), RangeSegment(0,N));

  auto knl1 = make_kernel<KPOL>(segment, l);
  auto knl2 = make_kernel<KPOL>(segment, l);

  auto d = create_layout_decisions<2>(tie(a));

  auto m0 = d.moments[0];
  auto m1 = d.moments[1];
  auto m2 = d.moments[2];

  auto la0 = get<0>(m0);
  auto la1 = get<0>(m1);
  auto la2 = get<0>(m2);

  auto p0 = layout_to_perm_array<2>(la0);
  auto p1 = layout_to_perm_array<2>(la1);
  auto p2 = layout_to_perm_array<2>(la2);

  EXPECT_EQ(perm, p0) << "Decision object starting layout matches original layout\n";
  EXPECT_EQ(perm, p1) << "Decision object layotu after first kernel matches original layout\n";
  EXPECT_EQ(perm, p2) << "Decision object layotu after second kernel matches original layout\n";
}

TEST(DecisionObjectTest, MultipleKernels) {
  using namespace RAJA;
  RAJA::idx_t N = 16;
  double * _a = new double[N*N];
  double * _b = new double[N*N];
 
  using KPOL = KernelPolicy<
    statement::For<0, loop_exec,
      statement::For<1, loop_exec,
        statement::Lambda<0>
      >
    >
  >;
  
  using VIEW2 = View<double, RAJA::Layout<2>>;
  
  std::array<camp::idx_t, 2> sizes{{N,N}};
  std::array<camp::idx_t, 2> perm{{1,0}};
  Layout<2> flipped = make_permuted_layout(sizes, perm);
  VIEW2 a(_a, flipped);
  VIEW2 b(_b, flipped);

  auto l = [&](auto i, auto j) {
    a(i,j) = 1;
  };

  auto segment = make_tuple(RangeSegment(0,N), RangeSegment(0,N));

  auto knl1 = make_kernel<KPOL>(segment, l);
  auto knl2 = make_kernel<KPOL>(segment, l);

  auto d2 = create_layout_decisions<2>(tie(a,b));
  auto d1 = create_layout_decisions<1>(tie(a,b));

  

}
TEST(DecisionObjectTest, FormatDecisions) {
  using namespace RAJA;
  RAJA::idx_t N = 16;
  double * _a = new double[N*N];
  double * _b = new double[N*N];
 
  using KPOL = KernelPolicy<
    statement::For<0, loop_exec,
      statement::For<1, loop_exec,
        statement::Lambda<0>
      >
    >
  >;
  
  using VIEW2 = View<double, RAJA::Layout<2>>;
  VIEW2 a(_a, N, N);
  VIEW2 b(_b, N, N);

  FormatDecisions<2, VIEW2&, VIEW2&> d(tie(a,b), camp::idx_seq<0,1>{});

  std::array<camp::idx_t, 2> sizes{{N,N}};
  std::array<camp::idx_t, 2> perm{{1,0}};
  Layout<2> flipped = make_permuted_layout(sizes, perm);

  d.change_to_before(a, flipped, 0);
  d.change_to_before(b, flipped, 1);

  EXPECT_EQ(d.user_decided[0][0], true);
  EXPECT_EQ(d.user_decided[1][1], true);
  EXPECT_EQ(d.user_decided[0][1], false);
  EXPECT_EQ(d.user_decided[1][0], false);
    
}
REGISTER_TYPED_TEST_SUITE_P(DecisionObjectTest, DecisionObject);
