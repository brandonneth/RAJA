#include "RAJA_test-base.hpp"

TEST(IslDecisionTest, DecisionVariableReps) {
  using namespace RAJA;
  std::string viewName = "name";
  auto perm = std::vector<idx_t>{{0,1}};
  auto d1 = new FormatDecision<2>("name", 0, perm);
  
  EXPECT_EQ(d1->rep(), "name_f0_0_1");
  
  auto c1 = new ConversionDecision<2>("name", 1, perm, perm);

  EXPECT_EQ(c1->rep(), "name_conv1_in_0_1_out_0_1");

  auto allvars = all_decision_variables<2>(viewName, 2);
  
  auto varlist = allvars->variables;

  std::vector<std::string> varlist_rep;
  for(auto v : varlist) {
    varlist_rep.push_back(v->rep());
  }

  std::vector<std::string> expected_reps;
  expected_reps.push_back("name_f0_0_1");
  expected_reps.push_back("name_f0_1_0");
  expected_reps.push_back("name_f1_0_1");
  expected_reps.push_back("name_f1_1_0");
  expected_reps.push_back("name_f2_0_1");
  expected_reps.push_back("name_f2_1_0");
  expected_reps.push_back("name_f3_0_1");
  expected_reps.push_back("name_f3_1_0");
  expected_reps.push_back("name_conv1_in_1_0_out_1_0");
  expected_reps.push_back("name_conv1_in_1_0_out_0_1");
  expected_reps.push_back("name_conv1_in_0_1_out_1_0");
  expected_reps.push_back("name_conv1_in_0_1_out_0_1");
  expected_reps.push_back("name_conv2_in_1_0_out_1_0");
  expected_reps.push_back("name_conv2_in_1_0_out_0_1");
  expected_reps.push_back("name_conv2_in_0_1_out_1_0");
  expected_reps.push_back("name_conv2_in_0_1_out_0_1");
  expected_reps.push_back("name_conv3_in_1_0_out_1_0");
  expected_reps.push_back("name_conv3_in_1_0_out_0_1");
  expected_reps.push_back("name_conv3_in_0_1_out_1_0");
  expected_reps.push_back("name_conv3_in_0_1_out_0_1");

  std::sort(varlist_rep.begin(), varlist_rep.end());
  std::sort(expected_reps.begin(), expected_reps.end());

  EXPECT_EQ(varlist_rep, expected_reps);
}

template <typename EXEC_POLICY>
void IslDecisionTestImpl() {

  return;
  camp::idx_t N = 1024;
  using namespace RAJA;
  
  View<double, Layout<2>> a(new double[N*N]);
  
  using KPOL = KernelPolicy<
    statement::For<0,loop_exec,
      statement::For<1, loop_exec,
        statement::Lambda<0>
      >
    >
  >;
  auto range = RangeSegment(0,N);
  auto segments = make_tuple(range, range);

  auto lambda = [&](auto i, auto j) { a(j,i) = 0; } ;
  auto knl = make_kernel<KPOL>(segments, lambda);

  auto accesses = knl.execute_symbolically();
  auto viewName = get_array_name(accesses[0]);

  std::cout << "View name from access: " << viewName << "\n";

  std::vector<camp::idx_t> perm1 = {0,1};
  std::vector<camp::idx_t> perm2 = {1,0};


  auto allvars = all_decision_variables<2>(viewName, 1);
  std::cout << "all vars: " << allvars->rep() << "\n";

  auto allconsts = all_constraints<2>(viewName, 1);
   
  std::cout << "all constraints: " << allconsts->rep() << "\n";

  

  isl_ctx * ctx = isl_ctx_alloc();
  isl_printer * p = isl_printer_to_file(ctx, stdout);
  p = isl_printer_set_suffix(p, "\n");

  isl_union_set * space = decision_space(ctx, allvars, allconsts);

  std::cout << "\nDecision Space\n";
  p = isl_printer_print_union_set(p,space);
  std::cout << "\n";
  Expression<2> * e = new MultExpression<2>(new Constant<2>(0), new Constant<2>(0));
  isl_union_map * obj = objective_function(ctx, allvars, allconsts, e);

  std::cout << "\nObjective Function Map\n";
  p = isl_printer_print_union_map(p, obj);

  std::cout << "\n";
  isl_union_set * result = isl_union_set_apply(space, obj);

  isl_union_set * best = isl_union_set_lexmin(result);
  isl_point * point = isl_union_set_sample_point(best);
  std::cout << "Selected Point\n";
  
  p = isl_printer_print_point(p,point);
  std::cout << "\n";
  

  for(int i = 0; i < allvars->variables.size(); i++) {
    std::cout << allvars->variables[i]->rep() << ": ";
    isl_val * v = isl_point_get_coordinate_val(point, isl_dim_all, i+1);
    p = isl_printer_print_val(p, v);
    std::cout << "\n";
  }
  std::cout << "\n";


  isl_ctx_free(ctx);


  
}


TEST(IslDecisionTest, MultipleKernelsOneView) {
  using namespace RAJA;
  using KPOL = KernelPolicy<
    statement::For<0, loop_exec,
      statement::For<1, loop_exec,
        statement::Lambda<0>
      >
    >
  >; 
  using VIEW2 = View<double, RAJA::Layout<2>>;

  RAJA::idx_t N = 16; 

  double * _a = new double[N*N];
  double * _b = new double[N*N];
  
  std::array<camp::idx_t, 2> sizes{{N,N}};
  std::array<camp::idx_t, 2> perm{{1,0}};
  Layout<2> flipped = make_permuted_layout(sizes, perm);
  Layout<2> normal = make_permuted_layout(sizes, {{0,1}});
  VIEW2 a(_a, flipped);
  VIEW2 b(_b, flipped);

  auto segs = make_tuple(RangeSegment(0,N), RangeSegment(0,N));

  
  auto l1 = [&](auto i, auto j) {
    a(i,j) = b(i,j) * 2;
  };
  auto l2 = [&](auto i, auto j) {
    a(i,j) = b(i,j) * 2;
  };
  auto l3 = [&](auto i, auto j) {
    a(i,j) = b(i,j) * 2;
  };

  auto knl1 = make_kernel<KPOL>(segs, l1);
  auto knl2 = make_kernel<KPOL>(segs, l2);
  auto knl3 = make_kernel<KPOL>(segs, l3);

  auto decisions1 = format_decisions(tie(a,b), knl1,knl2,knl3);

  decisions1.set_format_after(a, normal, knl1);
  std::cout << "view a: ";
  std::cout << decisions1.finalize_view(a) << "\n";
  std::cout << "view b: ";
  std::cout << decisions1.finalize_view(b) << "\n";
  
}



TYPED_TEST_SUITE_P(IslDecisionTest);
template <typename T>
class IslDecisionTest : public ::testing::Test
{
};

TYPED_TEST_P(IslDecisionTest, IslDecision)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  IslDecisionTestImpl<EXEC_POLICY>();

}




REGISTER_TYPED_TEST_SUITE_P(IslDecisionTest, IslDecision);
