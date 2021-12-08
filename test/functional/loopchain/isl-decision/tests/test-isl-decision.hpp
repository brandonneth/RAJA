#include "RAJA_test-base.hpp"



template <typename EXEC_POLICY>
void IslDecisionTestImpl() {

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
  isl_union_map * obj = objective_function(ctx, allvars, allconsts);

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
