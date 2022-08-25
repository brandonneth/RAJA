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
  
  auto constraints = base_constraints<2>("arr0", 1);
  
  
  isl_ctx * ctx = isl_ctx_alloc();
  isl_printer * pp = isl_printer_to_file(ctx, stdout);



  std::string aff_str = "{ [x,y] -> [100*x]}";
  isl_aff * aff = isl_aff_read_from_str(ctx, aff_str.c_str());
  std::cout << "\n\n";
  isl_printer_print_aff(pp, aff);
  std::cout << "\n\n";
  isl_qpolynomial * qpoly = isl_qpolynomial_from_aff(aff);
  isl_printer_print_qpolynomial(pp, qpoly);
  std::string poly_str = "{ [x, y] -> 100 * x * y }";
  isl_union_pw_qpolynomial * p = isl_union_pw_qpolynomial_read_from_str(ctx, poly_str.c_str());

  
  std::cout << "\n";
  isl_printer_print_union_pw_qpolynomial(pp, p);
  std::cout << "\n";
  


  isl_union_set * space = decision_space(ctx, modelVars, constraints);

  isl_printer_print_union_set(pp, space);
  std::cout << "\n\n";

  

  SumExpression<2> * objective_function = new SumExpression<2>();

  const std::vector<std::array<camp::idx_t, 2>> all_perms = enumerate_perms<2>();
  constexpr int numKnls = 1;
  //conversion component
  for(int i = 0; i < numKnls+1; i++) {
    for(auto in_perm : all_perms) {
      for(auto out_perm : all_perms) {
        auto in_var = new FormatDecision<2>("arr0", i, in_perm);
        auto out_var = new FormatDecision<2>("arr0", i+1, out_perm);
        auto estimated_cost = estimate_conversion_cost<2>(in_perm, out_perm);
        auto c = new Constant<2>(estimated_cost);
        MultExpression<2> * m = new MultExpression<2>(in_var, out_var);
        MultExpression<2> * term = new MultExpression<2>(m, c);
        objective_function->add(term);
      }
    }
  }
  std::cout << "conversion costs: ";
  std::cout << objective_function->rep() << "\n";
  //use component


  auto estimated_cost = choices1.template estimate_use_costs<0>(v1);

  objective_function->add(estimated_cost);
  


  std::cout << "with use costs: " << objective_function->rep() << "\n\n";


  std::stringstream polynomial;
  polynomial << "{" << modelVars->rep() << " -> " << objective_function->rep() << " } ";

  std::cout << "Polynomial representation of objective function: " << polynomial.str() << "\n";
  isl_union_pw_qpolynomial * obj = isl_union_pw_qpolynomial_read_from_str(ctx, polynomial.str().c_str());

  struct OptStruct {
    isl_union_pw_qpolynomial * obj;
    isl_point * lowest_point;
    camp::idx_t lowest_val;
  } fn_struct;

  fn_struct.obj = obj;
  fn_struct.lowest_point = NULL;
  fn_struct.lowest_val = INT_MAX;

  auto fn = [](isl_point * pnt, void * user) {
    OptStruct * fn_struct = (OptStruct*) user;
    isl_union_pw_qpolynomial * obj = isl_union_pw_qpolynomial_copy(fn_struct->obj);
    isl_point * lowest_point = fn_struct->lowest_point;
    camp::idx_t lowest_val = fn_struct->lowest_val;

    std::cout << "Applying polynomial.\n";
    isl_val * val = isl_union_pw_qpolynomial_eval(obj, isl_point_copy(pnt));
    std::cout << "Applied polynomial.\n";
    camp::idx_t v =  isl_val_get_num_si(val);
    std::cout << "val:" << v << "\n";
    if ( v < lowest_val) {
      fn_struct->lowest_val = v;
      fn_struct->lowest_point = pnt;
    }
    return isl_stat_ok;
  };

  isl_union_set_foreach_point(space, fn, &fn_struct);

  std::cout << "Lowest val: " << fn_struct.lowest_val << "\n";
  std::cout << "Lowest point: ";
  isl_printer_print_point(pp, fn_struct.lowest_point);
  std::cout << "\n";
  
  
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
