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



  auto eq = one_format_constraint<2>(viewName, 0);
  
  std::cout << "format constraint: " << eq->rep() << "\n";

  auto conv = one_conversion_constraint<2>(viewName, 1);
  std::cout << "conv constraint: " << conv->rep() << "\n";

  auto inputs = conversion_input_constraint<2>(viewName, 0);

  std::cout << "conversion input constraint: " << inputs->rep() << "\n";

  auto outputs = conversion_output_constraint<2>(viewName, 1);

  std::cout << "conversion output constraint: " << outputs->rep() << "\n";

  auto binary = binary_choice_constraints<2>(viewName, 0);

  std::cout << "binary choice for idx=0: " << binary->rep() << "\n";

  auto binary1 = binary_choice_constraints<2>(viewName, 1);

  std::cout << "binary choice for idx=1: " << binary1->rep() << "\n";

  auto vars0 = decision_variables<2>(viewName, 0);
  auto vars1 = decision_variables<2>(viewName, 1);

  std::cout << "lists of decision variables: " << vars0->rep() << "\n" << vars1->rep() << "\n";

  auto allvars = all_decision_variables<2>(viewName, 1);
  std::cout << "all vars: " << allvars->rep() << "\n";

  auto allconsts = all_constraints<2>(viewName, 1);
   
  std::cout << "all constraints: " << allconsts->rep() << "\n";


  isl_ctx * ctx = isl_ctx_alloc();

  isl_printer * p = isl_printer_to_file(ctx, stdout);

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

  std::cout << "Selected Point\n";
  p = isl_printer_print_union_set(p,best);

  for(int i = 0; i < allvars->variables.size(); i++) {
    std::cout << i << ": " << allvars->variables[i]->rep() << "\n";
  }
  std::cout << "\n";


  isl_ctx_free(ctx);

/*
  auto conv_var = new ConversionDecision<2>(viewName, 1, perm1, perm2);

  std::cout << "Conversion Variable for which I will estimate the cost: " << conv_var->rep() << "\n";

  auto cost = estimate_cost(conv_var);
  
  std::cout << "Estimated cost: " << cost << "\n";
  
  auto form_var = new FormatDecision<2>(viewName, 1, perm1);

   std::cout << "Format Variable for which I will estimate the cost: " << form_var->rep() << "\n";

  auto fcost = estimate_cost(form_var);
  
  std::cout << "Estimated cost: " << fcost << "\n";

  

  /*
  FormatDecision<camp::idx_t> format0_1(viewName, 0, perm1);
  FormatDecision<camp::idx_t> format0_2(viewName, 0, perm2);
  FormatDecision<camp::idx_t> format1_1(viewName, 1, perm1);
  FormatDecision<camp::idx_t> format1_2(viewName, 1, perm2);
  FormatDecision<camp::idx_t> format2_1(viewName, 2, perm1);
  FormatDecision<camp::idx_t> format2_2(viewName, 2, perm2);

  std::cout << "Format decision variable names:\n" 
    << format0_1.name() << "\n"
    << format0_2.name() << "\n"
    << format1_1.name() << "\n"
    << format1_2.name() << "\n"
    << format2_1.name() << "\n"
    << format2_2.name() << "\n";


  ConversionDecision<camp::idx_t> conv1_1(viewName, 1, perm1, perm1);
  ConversionDecision<camp::idx_t> conv1_2(viewName, 1, perm1, perm2);
  ConversionDecision<camp::idx_t> conv1_3(viewName, 1, perm2, perm1);
  ConversionDecision<camp::idx_t> conv1_4(viewName, 1, perm2, perm2);

  ConversionDecision<camp::idx_t> conv2_1(viewName, 2, perm1, perm1);
  ConversionDecision<camp::idx_t> conv2_2(viewName, 2, perm1, perm2);
  ConversionDecision<camp::idx_t> conv2_3(viewName, 2, perm2, perm1);
  ConversionDecision<camp::idx_t> conv2_4(viewName, 2, perm2, perm2);


  std::cout << "Conversion decision variable names:\n" 
    << conv1_1.name() << "\n"
    << conv1_2.name() << "\n"
    << conv1_3.name() << "\n"
    << conv1_4.name() << "\n"
    << conv2_1.name() << "\n"
    << conv2_2.name() << "\n"
    << conv2_3.name() << "\n"
    << conv2_4.name() << "\n";

  OneFormatConstraint<2> fc0(viewName, 0);
  OneFormatConstraint<2> fc1(viewName, 1);
  OneFormatConstraint<2> fc2(viewName, 2);

  std::cout << "Format Exclusivity Constraints:\n"
    << fc0.rep() << "\n"
    << fc1.rep() << "\n"
    << fc2.rep() << "\n";

 
  OneConversionConstraint<2> cc1(viewName,1);
  OneConversionConstraint<2> cc2(viewName,2);

  std::cout << "Conversion Exclusivity Constraints:\n"
    << cc1.rep() << "\n"
    << cc2.rep() << "\n";

  ConversionInputConstraint<2> cic1(viewName, 1);
  ConversionInputConstraint<2> cic1(viewName, 2);

  std::cout << "Conversion Input Constraints:\n"
   << cic1.rep() << "\n"
   << cic2.rep() << "\n";
  */
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
