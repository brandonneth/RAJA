#include "RAJA_test-base.hpp"

#include <unordered_set>
template <typename EXEC_POLICY>
void ILPModelTestImpl() {
  
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

  auto choices1 = format_decisions(tie(v1), knl1);

  //testing representation for various components
  auto formatVariable = new FormatDecision<2>("view", 0, perm_01);
  EXPECT_EQ(formatVariable->rep(), "view_comp0_0_1");

  auto termList = new TermList<2>();
  termList->add(formatVariable);
  EXPECT_EQ(termList->rep(), "[ view_comp0_0_1 ]");
  
  termList->add(formatVariable);
  EXPECT_EQ(termList->rep(), "[ view_comp0_0_1 , view_comp0_0_1 ]");

  //testing that the model variables are what is expected
  auto modelVars = all_decision_variables<2>("arr0", 1);
  
  auto expectedVars = std::unordered_set<std::string>();
  std::vector<std::string> moments = {"comp0", "comp1", "comp2"};
  std::vector<std::string> perms = {"0_1", "1_0"};

  for(auto m : moments) {
    for( auto p : perms) {
      std::string v = "arr0_";
      v += m;
      v += "_";
      v += p;
      expectedVars.insert(v);
    }
  }
  EXPECT_EQ(expectedVars.size(), modelVars->variables.size());
  for(auto v : modelVars->variables) {
    std::string r = v->rep();
    EXPECT_EQ(expectedVars.count(r), 1) << "failed for variable " << r;
  }
  
  //make sure that something gets chosen.

  using VIEW = View<double, Layout<2>>;

  idx_t ni = 1024;
  idx_t nj = 1024;
  idx_t nk = 1024;
  idx_t nl = 1024;
  idx_t nm = 1024;
  VIEW A(new double[ni*nk], ni,nk);
  VIEW B(new double[nj*nk], nk,nj);
  VIEW C(new double[nj*nm], nj,nm);
  VIEW D(new double[nm*nl], nm,nl);
  VIEW E(new double[ni*nj], ni,nj);
  VIEW F(new double[nj*nl], nj,nl);
  VIEW G(new double[ni*nl], ni,nl);
  
  auto lam_init_e = [&](auto i, auto j) {
    E(i,j) = 0;
  };
  auto lam_init_f = [&](auto i, auto j) {
    F(i,j) = 0;
  };
  auto lam_init_g = [&](auto i, auto j) {
    G(i,j) = 0;
  };

  auto lam_comp1 = [&](auto i, auto j, auto k) {
    E(i,j) += A(i,k) * B(k,j);
  };

  auto lam_comp2 = [&](auto i, auto j, auto k) {
    F(i,j) += C(i,k) * D(k,j);
  };
  auto lam_comp3 = [&](auto i, auto j, auto k) {
    G(i,j) += E(i,k) * F(k,j);
  };

  using namespace RAJA;
  using INIT_POL = KernelPolicy<
    statement::For<0, omp_parallel_for_exec,
      statement::For<1, loop_exec,
        statement::Lambda<0>
      >
    >
      >;
    
  using COMP_POL = KernelPolicy<
    statement::For<0, omp_parallel_for_exec,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
          statement::Lambda<0>
        >
      >
    >
  >;
  auto eseg = make_tuple(RangeSegment(0,ni), RangeSegment(0,nj));
  auto fseg = make_tuple(RangeSegment(0,nj), RangeSegment(0,nl));
  auto gseg = make_tuple(RangeSegment(0,ni), RangeSegment(0,nl));
  
  auto seg1 = make_tuple(RangeSegment(0,ni), RangeSegment(0,nj), RangeSegment(0,nk));
  auto seg2 = make_tuple(RangeSegment(0,nj), RangeSegment(0,nl), RangeSegment(0,nm));
  auto seg3 = make_tuple(RangeSegment(0,ni), RangeSegment(0,nl), RangeSegment(0,nj));


  auto init_e = make_kernel<INIT_POL>(eseg, lam_init_e);
  auto comp1 = make_kernel<COMP_POL>(seg1, lam_comp1);
  auto init_f = make_kernel<INIT_POL>(fseg, lam_init_f);
  auto comp2 = make_kernel<COMP_POL>(seg2, lam_comp2);
  auto init_g = make_kernel<INIT_POL>(gseg, lam_init_g);
  auto comp3 = make_kernel<COMP_POL>(seg3, lam_comp3);

  std::array<idx_t, 2> bsizes{{nk,nj}};
  auto blayout = make_permuted_layout(bsizes, {{1,0}});
  std::array<idx_t, 2> dsizes{{nm,nl}};
  auto dlayout = make_permuted_layout(dsizes, {{1,0}});
  std::array<idx_t, 2> fsizes{{nj,nl}};
  auto flayout01 = make_permuted_layout(fsizes, {{0,1}});
  auto flayout10 = make_permuted_layout(fsizes, {{1,0}});
  
  

  auto dec = format_decisions(tie(B, D, F), init_e, comp1, init_f, comp2, init_g, comp3);
  
  dec.set_format_for(B, blayout, comp1);
  dec.set_format_for(D, blayout, comp2);
  dec.set_format_for(F, flayout01, init_f);
  dec.set_format_for(F, flayout10, comp3);
  
  auto decision_variables = all_decision_variables<2>(dec.view_name(B), 6);

  SumExpression<2> * use_cost = dec.template estimate_use_costs<0>(B);
  SumExpression<2> * conversion_cost = dec.template estimate_conversion_costs(B);
  SumExpression<2> * cost_function = new SumExpression<2>();
  cost_function->add(use_cost);
  cost_function->add(conversion_cost);

  auto base_constraints_ = base_constraints<2>(dec.view_name(B), 6); 
  auto extra_constraints = dec.template user_constraints<0>(B);
  auto all_constraints = new AndConstraint<2>(base_constraints_, extra_constraints);

  isl_ctx * ctx2 = isl_ctx_alloc();

  isl_union_set * all_choices = decision_space(ctx2, decision_variables, all_constraints);
  EXPECT_FALSE(isl_union_set_is_empty(all_choices)) << "No elements in the choices set.\n";
  
  isl_union_pw_qpolynomial * polynomial = objective_function(ctx2, decision_variables, cost_function);
  isl_point * choice = find_best_choice(all_choices, polynomial);

  EXPECT_FALSE(isl_point_is_void(choice)) << "Objective function did not generate any point.\n";

}

TYPED_TEST_SUITE_P(ILPModelTest);
template <typename T>
class ILPModelTest : public ::testing::Test
{
};

TYPED_TEST_P(ILPModelTest, FormatDecisions)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  ILPModelTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(ILPModelTest, FormatDecisions);
