#include "RAJA_test-base.hpp"

#include <unordered_set>
template <typename EXEC_POLICY>
void NestingOrderTestImpl() {
  using namespace RAJA;
  using pol1 = KernelPolicy<
    statement::For<0,loop_exec,
      statement::For<1,loop_exec,
        statement::For<2,loop_exec,
          statement::Lambda<0>,
          statement::Lambda<1>
        >,
        statement::Lambda<2>
      >
    >
  >;
  std::vector<idx_t> correct1_0{{0,1,2}};
  std::vector<idx_t> correct1_1 = correct1_0;
  std::vector<idx_t> correct1_2{{0,1}};
  auto order1_0 = policy_to_nesting_order<pol1,0>(); 
  auto order1_1 = policy_to_nesting_order<pol1,1>(); 
  auto order1_2 = policy_to_nesting_order<pol1,2>(); 

  EXPECT_EQ(order1_0, correct1_0);
  EXPECT_EQ(order1_1, correct1_1);
  EXPECT_EQ(order1_1, correct1_1);

  using pol2 = KernelPolicy<
    statement::For<2,loop_exec,
      statement::For<0,loop_exec,
        statement::Lambda<0>
      >,
      statement::For<1,loop_exec,
        statement::Lambda<1>
      >,
      statement::For<0,loop_exec,
        statement::Lambda<2>
      >
    >
  >;
  std::vector<idx_t> correct2_0{{2,0}};
  std::vector<idx_t> correct2_1{{2,1}};
  std::vector<idx_t> correct2_2{{2,0}};
  auto order2_0 = policy_to_nesting_order<pol2,0>(); 
  auto order2_1 = policy_to_nesting_order<pol2,1>(); 
  auto order2_2 = policy_to_nesting_order<pol2,2>(); 

  EXPECT_EQ(order2_0, correct2_0);
  EXPECT_EQ(order2_1, correct2_1);
  EXPECT_EQ(order2_1, correct2_1);

}
template <typename EXEC_POLICY>
void ImperfectILPModelTestImpl() {
  using namespace RAJA;
  using MEANPOL = KernelPolicy<
        statement::For<0, omp_parallel_for_exec,
          statement::Lambda<0, Segs<0>>,
          statement::For<1, loop_exec,
            statement::Lambda<1, Segs<0,1>>
          >,
          statement::Lambda<2, Segs<0>>
        >
      >;
  std::array<RAJA::idx_t, 2> perm_01 {{0,1}};
  std::array<RAJA::idx_t, 2> perm_10 {{1,0}};
  std::array<RAJA::idx_t, 1> perm_0 {{0}};

  using VIEW2 = RAJA::View<double, RAJA::Layout<2,RAJA::idx_t>>;
  using VIEW1 = RAJA::View<double, RAJA::Layout<1,RAJA::idx_t>>;

  idx_t ni = 10;
  idx_t nj = 7;

  double * _mean = new double[ni];
  double * _data = new double[ni*nj];
  VIEW1 mean(_mean, make_permuted_layout({{nj}}, perm_0));
  VIEW1 stddev(new double[nj], make_permuted_layout({{nj}}, perm_0));
  VIEW2 data(_data, make_permuted_layout({{ni,nj}}, perm_01));
   RAJA::View<Real_type, RAJA::Layout<2>> corr(new double[nj*nj], nj, nj);\
  for(int i = 0; i < ni; i++) {
    for(int j = 0; j < nj; j++) {
      data(i,j) = j;
    }
  } 


      auto correlation_lam1 = [&](auto j) {mean(j) = 0.0;};

      auto correlation_lam2 = [&](auto j, auto i) { mean(j) += data(i,j);  };

      auto correlation_lam3 = [&](auto j) {mean(j) /= ni;};

      auto correlation_lam4 = [&](auto j) {stddev(j) = 0.0;};

      auto correlation_lam5 = [&](auto j, auto i) { stddev(j) += (data(i,j) - mean(j)) * (data(i,j) - mean(j));  };

      auto correlation_lam6 = [&](auto j) {
        stddev(j) /= ni;
        stddev(j) = std::sqrt(stddev(j));
        stddev(j) = stddev(j) <= 0.001 ? 1.0 : stddev(j);
      };

      auto correlation_lam7 = [&](auto i, auto j) {
        data(i,j) -= mean(j);
        data(i,j) /= ni * stddev(j);
      };

      auto correlation_lam8 = [&](auto i0) {
         corr(i0,i0) = 1.0;
      };

      auto correlation_lam9 = [&](auto i0, auto i1) {
         corr(i0,i1) = 0.0;
      };

      auto correlation_lam10 = [&](auto i, auto j, auto k) {
         corr(i,j) += data(k,i) * data(k,j);
      };

      auto correlation_lam11 = [&](auto i, auto j) {
        corr(j,i) = corr(i,j);
      };

      using MEANPOL = KernelPolicy<
        statement::For<0, omp_parallel_for_exec,
          statement::Lambda<0, Segs<0>>,
          statement::For<1, loop_exec,
            statement::Lambda<1, Segs<0,1>>
          >,
          statement::Lambda<2, Segs<0>>
        >
      >;
      auto mean_segs = make_tuple(RangeSegment(0,nj), RangeSegment(0,ni));
      auto mean_knl = make_kernel<MEANPOL>(mean_segs, correlation_lam1, correlation_lam2, correlation_lam3);
      using DEVPOL=MEANPOL;
      auto dev_segs = make_tuple(RangeSegment(0,nj), RangeSegment(0,ni));
      auto dev_knl = make_kernel<DEVPOL>(dev_segs, correlation_lam4, correlation_lam5, correlation_lam6);

      using CENTERPOL = KernelPolicy<
        statement::For<0,omp_parallel_for_exec,
          statement::For<1,loop_exec,
            statement::Lambda<0>
          >
        >
      >;
      auto centersegs = make_tuple(RangeSegment(0,ni), RangeSegment(0,nj));
      auto center_knl = make_kernel<CENTERPOL>(centersegs, correlation_lam7);
       using CALCPOL = KernelPolicy<
        statement::For<0,loop_exec,
          statement::Lambda<0,Segs<0>>,
          statement::For<1, loop_exec,
            statement::Lambda<1,Segs<0,1>>,
            statement::For<2,loop_exec,
              statement::Lambda<2,Segs<0,1,2>>
            >,
            statement::Lambda<3,Segs<0,1>>
          >
        >
      >;
      auto calc_seg1 = make_symbolic_segment(0,nj-1);
      auto calc_seg2 = make_symbolic_segment(calc_seg1 + 1, nj);
      auto calc_seg3 = make_symbolic_segment(0,ni);
      auto calc_segs = make_tuple(calc_seg1, calc_seg2, calc_seg3);
      auto calc_knl = make_kernel<CALCPOL>(calc_segs, correlation_lam8, correlation_lam9, correlation_lam10, correlation_lam11);



     auto dec = format_decisions(tie(data), mean_knl, dev_knl, center_knl, calc_knl);
      auto comp = dec.finalize();

   std::cout << "finalized computation for correlation: ";
   dec.print_to_stream(std::cout);
   std::cout << "\n";
   comp();

  auto dec2 = format_decisions_ref(tie(data), tie(mean_knl, dev_knl, center_knl, calc_knl));
  dec2.set_format_for(data, {1,0}, mean_knl);
  dec2.set_format_for(data, {1,0}, calc_knl);
  auto comp2 = dec2.finalize();
  std::cout << "\nfor reference using one: ";
  dec2.print_to_stream(std::cout);
  std::cout << "\n";
  comp2();
}

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
  std::cout << "made choices1\n";
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
  std::cout << "testing modelVars and expectedVars\n";
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
  dec.set_format_for(D, dlayout, comp2);
  dec.set_format_for(F, flayout01, init_f);
  dec.set_format_for(F, flayout10, comp3);
  std::cout << "set the variables for the 3mm situation\n"; 
  auto decision_variables = all_decision_variables<2>(dec.view_name(B), 6);
  std::cout << "created the decision variable list\n";
  SumExpression<2> * use_cost = dec.template estimate_use_costs<0>(B);
  std::cout << "created use costs\n" << std::flush;
  SumExpression<2> * conversion_cost = dec.template estimate_conversion_costs(B);
  std::cout << "created use and conversion costs\n";
  SumExpression<2> * cost_function = new SumExpression<2>();
  cost_function->add(use_cost);
  cost_function->add(conversion_cost);
  std::cout << "created and added to the cost function\n";
  auto base_constraints_ = base_constraints<2>(dec.view_name(B), 6); 
  auto extra_constraints = dec.template user_constraints<0>(B);
  auto all_constraints = new AndConstraint<2>(base_constraints_, extra_constraints);

  isl_ctx * ctx2 = isl_ctx_alloc();

  isl_union_set * all_choices = decision_space(ctx2, decision_variables, all_constraints);
  std::cout << "created the decision space\n";
  EXPECT_FALSE(isl_union_set_is_empty(all_choices)) << "No elements in the choices set.\n";
  std::cout << "tested all_choices\n";
  isl_union_pw_qpolynomial * polynomial = objective_function(ctx2, decision_variables, cost_function);
  isl_point * choice = dec.find_best_choice(all_choices, polynomial);
  std::cout << "tested polynomial\n";
  EXPECT_FALSE(isl_point_is_void(choice)) << "Objective function did not generate any point.\n";
  std::cout << "finished the first ilpmodel test function\n";
}

TYPED_TEST_SUITE_P(ILPModelTest);
template <typename T>
class ILPModelTest : public ::testing::Test
{
};

TYPED_TEST_P(ILPModelTest, FormatDecisions)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  NestingOrderTestImpl<EXEC_POLICY>();
  ILPModelTestImpl<EXEC_POLICY>();
  
  ImperfectILPModelTestImpl<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(ILPModelTest, FormatDecisions);
