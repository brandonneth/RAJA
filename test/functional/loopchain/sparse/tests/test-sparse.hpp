#include "RAJA_test-base.hpp"
#include <string>
template <typename EXEC_POLICY>
void SparseTestImpl() {
  
  std::vector<idx_t> r{0,1,1,2,2,2,2,3,5,5};
  std::vector<idx_t> c{3,1,4,0,2,3,5,5,2,4};
  std::vector<double>v{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};
  SparseView view1(r,c,v);

  for(int i = 0; i < r.size(); i++ ) {
    EXPECT_EQ(view1(r(i), c(i)), v(i)) << "Failed for i = " << i << ".";
  }

}

TYPED_TEST_SUITE_P(SparseTest);
template <typename T>
class SparseTest : public ::testing::Test
{
};

TYPED_TEST_P(SparseTest, SymEval)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  SparseTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(SparseTest, SymEval);
