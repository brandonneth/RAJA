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
    statement::For<0,loop_exec,
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
  EXPECT_TRUE(choices.isl_time > 0);
  EXPECT_TRUE(choices.setup_time > 0);
  EXPECT_TRUE(choices.solve_time > 0);

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
