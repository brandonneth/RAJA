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
