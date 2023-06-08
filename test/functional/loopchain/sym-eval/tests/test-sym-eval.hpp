#include "RAJA_test-base.hpp"
#include <string>
template <typename EXEC_POLICY>
void SymEvalTestImpl() {
  
  double * _a = new double[100];
  double * _b = new double[100];

  using VIEW_TYPE = RAJA::View<double, RAJA::Layout<1, int, 1>>;

  VIEW_TYPE a(_a, 100);
  VIEW_TYPE b(_b, 100);

  for(int i = 0; i < 100; i++) {
    _a[i] = 10;
  }

  //reads, writes, and array matching
  auto forall1 = RAJA::make_forall<RAJA::seq_exec>(RAJA::RangeSegment(0,100), [&](auto i) {b(i) = a(i);});
  auto symbolicAccesses = forall1.execute_symbolically();
  EXPECT_EQ(2, symbolicAccesses.size());

  for(auto acc : symbolicAccesses) {
    if(acc.view == _a) {
      EXPECT_TRUE(acc.isRead) << "sym eval did not gather read of a\n";
      EXPECT_TRUE(!acc.isWrite) << "sym eval gathered extraneous write to a\n";
    } else if (acc.view == _b) {
      EXPECT_TRUE(!acc.isRead) << "sym eval gathere extraneous read from b\n";
      EXPECT_TRUE(acc.isWrite) << "sym eval did not gather write to b\n";
    } else{ 
      FAIL() << "Access view is not a or b\n";
    }
  }
  
  //access expression
  
  auto forall2 = RAJA::make_forall<RAJA::seq_exec>(RAJA::RangeSegment(0,100), [&](auto i) {a(i+2) = 0;});

  auto sym = forall2.execute_symbolically();

  EXPECT_EQ(1, sym.size());

  auto acc = sym[0];

  std::string expected_access_string = "i0+2";
  
  EXPECT_TRUE(acc.access_string().find(expected_access_string) != std::string::npos) << "Access string for a(i+2) was unexpected. It was: " << acc.access_string() << ", which doesn't contain " << expected_access_string << "\n";

  //TODO: 2d stuff
  //TODO: layout-to-perm
  free(_a);
  free(_b);
}



template <typename EXEC_POLICY>
void SymEvalConditionalTestImpl() {
  const size_t N = 10;
  double * _a = new double[N*N];
  double * _b = new double[N*N];

  using VIEW_TYPE = RAJA::View<double, RAJA::Layout<2>>;
  using namespace RAJA; 

  VIEW_TYPE a(_a, N, N);
  VIEW_TYPE b(_b, N, N);
  
  auto equals = [&](auto i, auto j) {
    if (i == j) {
      a(i,j) = b(i,j);
    }
  };

  auto not_equals = [&](auto i, auto j) {
    if (i != j) {
      a(i,j) = b(i,j);
    }
  };

  auto segs = make_tuple(RangeSegment(0,N),RangeSegment(0,N));
  using namespace statement;
  using pol = KernelPolicy<For<0,loop_exec,For<1,loop_exec,Lambda<0>>>>;
  auto eq_knl = make_kernel<pol>(segs, equals);
  auto neq_knl = make_kernel<pol>(segs, not_equals);
  
  auto eq_accesses01 = eq_knl.execute_symbolically(0,1);
  auto eq_accesses00 = eq_knl.execute_symbolically(0,0);
  
  EXPECT_EQ(eq_accesses01.size(), 0);
  EXPECT_EQ(eq_accesses00.size(), 2);
  
  auto neq_accesses01 = neq_knl.execute_symbolically(0,1);
  auto neq_accesses00 = neq_knl.execute_symbolically(0,0);
  
  EXPECT_EQ(neq_accesses01.size(), 2);
  EXPECT_EQ(neq_accesses00.size(), 0);
}

template <typename EP>
void SymEvalSparseTestImpl() {
  using namespace RAJA;

  auto N = 10;
  auto density = 0.01;
  auto A = make_random_sparse_view2<double>(N,density);
double temp = 0.0;

using POLICY = KernelPolicy<
  statement::For<0,loop_exec,
    statement::Lambda<0, Segs<0>>,
    statement::For<1,loop_exec,
      statement::Lambda<1>
    >,
    statement::Lambda<2,Segs<0>>
  >
>;
std::stringstream out;
auto lam1 = [&](auto i) {
  out << "1" << " " << i << "\n";
};
auto lam2 = [&](auto i, auto j) {
  if (j != i) {
    out << 2 << " " << i << " " << j << "\n";
    temp = A(i,j);
  }
};
auto lam3 = [&](auto i) {
  out << 3 << " " << i << "\n";
};

auto seg1 = RangeSegment(0,N);
auto seg2 = RangeSegment(0,N);
auto segs = make_tuple(seg1, seg2);

auto knl01 = make_sparse_kernel_sym<POLICY, 1>(make_tuple(0,1), segs, A, lam1, lam2, lam3);
auto knl00 = make_sparse_kernel_sym<POLICY, 1>(make_tuple(0,0), segs, A, lam1, lam2, lam3);

out.clear();
out.str("");

knl00();
std::string out00 = out.str();
std::cout << "out00\n" << out00;;


out.clear();
out.str("");
knl01();
std::string out01 = out.str();
out.clear();
out.str("");

EXPECT_NE(out01, out00);




}




TYPED_TEST_SUITE_P(SymEvalTest);
template <typename T>
class SymEvalTest : public ::testing::Test
{
};

TYPED_TEST_P(SymEvalTest, SymEval)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  std::cout << "Running symeval test";
  SymEvalTestImpl<EXEC_POLICY>();
  SymEvalConditionalTestImpl<EXEC_POLICY>();
  SymEvalSparseTestImpl<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SymEvalTest, SymEval);
