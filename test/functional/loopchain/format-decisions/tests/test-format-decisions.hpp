#include "RAJA_test-base.hpp"

template <typename EXEC_POLICY>
void FormatDecisionsTestImpl() {
  
  //permute_view
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


  auto data1 = new double[8*8];
  auto data2 = new double[8*8];
  auto l01 = RAJA::make_permuted_layout( {{8,8}}, perm_01 );
  auto l10 = RAJA::make_permuted_layout( {{8,8}}, perm_10 );

  VIEW2 v1(data1, l01);
  VIEW2 v2(data2, l10);


  auto segs = make_tuple(RangeSegment(0,8), RangeSegment(0,8));
  using KPOL = KernelPolicy<
    statement::For<0,EXEC_POLICY,
      statement::For<1,loop_exec,
        statement::Lambda<0>
      >
    >
  >;
  auto knl1 = make_kernel<KPOL>(segs, [&](auto i, auto j) {v1(i,j) = 2*i+j;});
  auto knl2 = make_kernel<KPOL>(segs, [&](auto i, auto j) {v2(i,j) = v1(i,j) * 2;});

  auto choices = format_decisions(camp::tie(v1,v2), knl1, knl2);
  
  //view and kernel index functions
  EXPECT_EQ(choices.view_index(v1), 0);
  EXPECT_EQ(choices.view_index(v2), 1);
  EXPECT_EQ(choices.kernel_index(knl1), 0);
  EXPECT_EQ(choices.kernel_index(knl2), 1);

  auto vector_01 = std::vector<RAJA::idx_t>{{0,1}};
  auto vector_10 = std::vector<RAJA::idx_t>{{1,0}};

  //get_layout_choice
  auto choice1 = choices.template get_layout_choice<0>(v1);
  EXPECT_EQ(choice1, vector_01);
  auto choice2 = choices.template get_layout_choice<0>(v2);
  EXPECT_EQ(choice2, vector_10);

  //set_format_for
  choices.set_format_for(v1, l10, knl1);
  choices.set_format_for(v1, l01, knl2);
  choices.set_format_for(v2, l01, knl1);
  choices.set_format_for(v2, l01, knl2);
  auto choice3 = choices.template get_layout_choice<1>(v1);
  EXPECT_EQ(choice3, vector_10);
  auto choice4 = choices.template get_layout_choice<2>(v1);
  EXPECT_EQ(choice4, vector_01);
  auto choice5 = choices.template get_layout_choice<1>(v2);
  EXPECT_EQ(choice5, vector_01);
  auto choice6 = choices.template get_layout_choice<2>(v2);
  EXPECT_EQ(choice6, vector_01);

  //set_output_format
  choices.set_output_format(v1, l10);
  choices.set_output_format(v2, l01);
  auto choice7 = choices.template get_layout_choice<3>(v1);
  EXPECT_EQ(choice7, vector_10);
  auto choice8 = choices.template get_layout_choice<3>(v2);
  EXPECT_EQ(choice8, vector_01);

  choices.lock();

  auto computation = choices.finalize();

  computation();

  EXPECT_EQ(layout_to_perm(v1.get_layout()), vector_10);
  EXPECT_EQ(layout_to_perm(v2.get_layout()), vector_01);
  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      EXPECT_EQ(v1(i,j), 2*i+j);
    }
  }
  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      EXPECT_EQ(v1(i,j) * 2, v2(i,j));
    }
  }

  auto times = choices.time_execution();
  EXPECT_TRUE(camp::get<0>(times) > 0);
  EXPECT_TRUE(camp::get<1>(times) > 0);
  
  //test simpler permutation description.
  auto choices2 = format_decisions(camp::tie(v1,v2), knl1, knl2);
  choices2.set_format_for(v1, {{1,0}}, knl1);
  choices2.set_format_for(v1, {{0,1}}, knl2);
  choices2.set_format_for(v2, {{0,1}}, knl1);
  choices2.set_format_for(v2, {{0,1}}, knl2);
  EXPECT_EQ(choices2.template get_layout_choice<1>(v1), vector_10);
  EXPECT_EQ(choices2.template get_layout_choice<2>(v1), vector_01);
  EXPECT_EQ(choices2.template get_layout_choice<1>(v2), vector_01);
  EXPECT_EQ(choices2.template get_layout_choice<2>(v2), vector_01);

  //set_output_format
  choices2.set_output_format(v1, {1,0});
  choices2.set_output_format(v2, {0,1});
  EXPECT_EQ(choices2.template get_layout_choice<3>(v1), vector_10);
  EXPECT_EQ(choices2.template get_layout_choice<3>(v2), vector_01);


  auto choices3 = format_decisions(camp::tie(v1,v2), knl1, knl2);
  choices3.set_format_for(v1, {{1,0}}, knl1);
  choices3.set_format_for(v1, {{0,1}}, knl2);
  choices3.set_format_for(v2, {{0,1}}, knl1);
  choices3.set_format_for(v2, {{0,1}}, knl2);
  EXPECT_EQ(choices3.template get_layout_choice<1>(v1), vector_10);
  EXPECT_EQ(choices3.template get_layout_choice<2>(v1), vector_01);
  EXPECT_EQ(choices3.template get_layout_choice<1>(v2), vector_01);
  EXPECT_EQ(choices3.template get_layout_choice<2>(v2), vector_01);

  //set_output_format
  choices3.set_output_format(v1, {1,0});
  choices3.set_output_format(v2, {0,1});
  EXPECT_EQ(choices3.template get_layout_choice<3>(v1), vector_10);
  EXPECT_EQ(choices3.template get_layout_choice<3>(v2), vector_01);
}


template <typename EXEC_POLICY>
void InliningTest() {
  using namespace RAJA;
  idx_t n = 100;
  using VIEW = View<double, Layout<2>>;
  VIEW A(new double[n*n], n,n);
  VIEW B(new double[n*n], n,n);
  VIEW C(new double[n*n], n,n);
  VIEW T(new double[n*n], n,n);
  VIEW D(new double[n*n], n,n);

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      T(i,j) = 0;
      D(i,j) = 0;
      A(i,j) = rand();
      B(i,j) = rand();
      C(i,j) = rand();
    }
  }

  using COMP_POL = KernelPolicy<
    statement::For<0, EXEC_POLICY,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
            statement::Lambda<0>
        >
      >
    >
  >;

  auto lam1 = [&](auto i, auto j, auto k) {
    T(i,j) += A(i,k) * B(k,j);
  };

  auto lam2 = [&](auto i, auto j, auto k) {
    D(i,j) += T(i,k) * C(k,j);
  };
  auto comp_seg = make_tuple(RangeSegment(0,n), RangeSegment(0,n), RangeSegment(0,n));

  auto comp1 = make_kernel<COMP_POL>(comp_seg, lam1);
  auto comp2 = make_kernel<COMP_POL>(comp_seg, lam2);

  auto fd = format_decisions(tie(B,C), comp1, comp2);
  fd.set_format_for(B,{1,0},comp1);
  fd.set_format_for(B,{1,0},comp2);
  fd.set_format_for(C,{0,1},comp1);
  fd.set_format_for(C,{1,0},comp2);
  fd.set_output_format(B,{1,0});
  fd.set_output_format(C,{1,0});
  fd.lock();

  auto times = fd.time_execution();
}


template <typename EPol>
void GetSeedPermsTest() {
  using namespace RAJA;
  idx_t n = 100;
  auto perm01 = std::array<idx_t, 2>{0,1};
  auto perm10 = std::array<idx_t, 2>{1,0};
  auto perm012 = std::array<idx_t, 3>{0,1,2};
  auto perm210 = std::array<idx_t, 3>{2,1,0};

  View<double, Layout<2>> A(new double[n*n], n,n);
  View<double, Layout<3>> B(new double[n*n*n], n,n,n);
  View<double, Layout<2>> C(new double[n*n], make_permuted_layout(std::array<idx_t, 2>{n,n}, perm10));

  auto knl1 = make_forall<EPol>(RangeSegment(0,n),[&](auto i) {B(i,i,i) = A(i,i);});
  auto knl2 = make_forall<EPol>(RangeSegment(0,n),[&](auto i) {B(i+1,i,i) = C(i,i);});

  auto dec = format_decisions(tie(A,B,C), knl1, knl2);

  //at the start it should be just the data's original layout.
  auto permsA1 = dec.get_seed_perms(A);
  EXPECT_EQ(1, permsA1.size());
  for(auto perm: permsA1) {
    auto expected = std::set<std::array<idx_t, 2>>{perm01};
    EXPECT_EQ(1, expected.count(perm)); 
  }
  auto permsB1 = dec.get_seed_perms(B);
  EXPECT_EQ(1, permsB1.size());
  for(auto perm: permsB1) {
    auto expected = std::set<std::array<idx_t, 3>>{perm012};
    EXPECT_EQ(1, expected.count(perm)); 
  }
  auto permsC1 = dec.get_seed_perms(C);
  EXPECT_EQ(1, permsC1.size());
  for(auto perm: permsC1) {
    auto expected = std::set<std::array<idx_t, 2>>{perm10};
    EXPECT_EQ(1, expected.count(perm)); 
  }

  dec.set_format_for(A, {1,0}, knl1);
  dec.set_format_for(B, {2,1,0}, knl2);
  dec.set_format_for(C, {1,0}, knl2);

  //it should include all formats that are registered.
  auto permsA2 = dec.get_seed_perms(A);
  EXPECT_EQ(2, permsA2.size());
  for(auto perm: permsA2) {
    auto expected = std::set<std::array<idx_t, 2>>{perm01, perm10};
    EXPECT_EQ(1, expected.count(perm)); 
  }
  auto permsB2 = dec.get_seed_perms(B);
  EXPECT_EQ(2, permsB2.size());
  for(auto perm: permsB2) {
    auto expected = std::set<std::array<idx_t, 3>>{perm012, perm210};
    EXPECT_EQ(1, expected.count(perm)); 
  }
  auto permsC2 = dec.get_seed_perms(C);
  EXPECT_EQ(1, permsC2.size());
  for(auto perm: permsC2) {
    auto expected = std::set<std::array<idx_t, 2>>{perm10};
    EXPECT_EQ(1, expected.count(perm)); 
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

  FormatDecisionsTestImpl<RAJA::loop_exec>();
  FormatDecisionsTestImpl<RAJA::omp_parallel_for_exec>();
  InliningTest<RAJA::omp_parallel_for_exec>();

  GetSeedPermsTest<RAJA::loop_exec>();
}

REGISTER_TYPED_TEST_SUITE_P(FormatDecisionsTest, FormatDecisions);
