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


TEST(SymEvalTest, NormalizeAccess) {
  using namespace RAJA;
  using VIEW = View<double, Layout<3>>;

  int N = 8;
  auto segs = make_tuple(RangeSegment(0,N), RangeSegment(0,N), RangeSegment(0,N));
  VIEW a(new double[N*N*N], N,N,N);

  Layout<3> layout1 = make_permuted_layout<3>({{N,N,N}}, {{0,1,2}});
  Layout<3> layout2 = make_permuted_layout<3>({{N,N,N}}, {{1,0,2}});

  using KPOL1 = KernelPolicy<
    statement::For<0,loop_exec,
      statement::For<1,loop_exec,
        statement::For<2,loop_exec,
          statement::Lambda<0>
        >
      >
    >
  >;
  using KPOL2 = KernelPolicy<
    statement::For<2,loop_exec,
      statement::For<0,loop_exec,
        statement::For<1,loop_exec,
          statement::Lambda<0>
        >
      >
    >
  >;

  
  a.set_layout(layout1);
  auto knl11 = make_kernel<KPOL1>(segs, [=](auto i, auto j, auto k) {a(i,j,k) = 1;});
  auto access11 = knl11.execute_symbolically().at(0);
  auto n11 = knl11.normalize_access(access11);
  EXPECT_EQ(n11[0], 0);
  EXPECT_EQ(n11[1], 1);
  EXPECT_EQ(n11[2], 2);

  a.set_layout(layout2);
  auto knl12 = make_kernel<KPOL1>(segs, [=](auto i, auto j, auto k) {a(i,j,k) = 1;});
  auto access12 = knl12.execute_symbolically().at(0);
  auto n12 = knl12.normalize_access(access12);
  EXPECT_EQ(n12[0], 1);
  EXPECT_EQ(n12[1], 0);
  EXPECT_EQ(n12[2], 2);

  a.set_layout(layout1);
  auto knl21 = make_kernel<KPOL2>(segs, [=](auto i, auto j, auto k) {a(i,j,k) = 1;});
  auto access21 = knl21.execute_symbolically().at(0);
  auto n21 = knl21.normalize_access(access21);
  EXPECT_EQ(n21[0], 1);
  EXPECT_EQ(n21[1], 2);
  EXPECT_EQ(n21[2], 0);


  a.set_layout(layout2);
  auto knl22 = make_kernel<KPOL2>(segs, [=](auto i, auto j, auto k) {a(i,j,k) = 1;});
  auto access22 = knl22.execute_symbolically().at(0);
  auto n22 = knl22.normalize_access(access22);
  EXPECT_EQ(n22[0], 2);
  EXPECT_EQ(n22[1], 1);
  EXPECT_EQ(n22[2], 0);

}


TYPED_TEST_SUITE_P(SymEvalTest);
template <typename T>
class SymEvalTest : public ::testing::Test
{
};

TYPED_TEST_P(SymEvalTest, SymEval)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  SymEvalTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(SymEvalTest, SymEval);
