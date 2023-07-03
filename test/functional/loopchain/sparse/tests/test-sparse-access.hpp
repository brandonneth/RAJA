#include "RAJA_test-base.hpp"
#include <string>
#include <vector>





template <typename EXEC_POLICY>
void SparseAccessTestImpl() {
   using namespace RAJA;
  std::vector<camp::idx_t> r{0,1,1,2,2,2,2,3,5,5};
  std::vector<camp::idx_t> c{3,1,4,0,2,3,5,5,2,4};
  std::vector<double>v{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};


  //SparseView<2,double> view1 = make_sparse_view<double>(r,c,v);
  auto view1 = make_sparse_view<double,2,COO>(r,c,v);

  for(size_t i = 0; i < r.size(); i++ ) {
    EXPECT_EQ(view1(r[i], c[i]), v[i]) << "Failed for i = " << i << ".";
  }



  using POL = KernelPolicy<
    statement::For<0, loop_exec,
      statement::For<1, loop_exec,
        statement::Lambda<0>
      >
    >
  >; //perfect_nest<0,1>::policy;

  auto seg = RangeSegment(0,6);
  auto segs = make_tuple(seg, seg);

  double sum = 0.0;
  auto knl = make_kernel<POL>(segs, [&](auto i, auto j) {
    sum += view1(i,j);
  });

  auto accesses = knl.execute_symbolically();

  EXPECT_EQ(accesses.size(), 1);

  auto access = accesses[0];

  EXPECT_EQ(access.view, view1.get_data());
  EXPECT_EQ(true, access.is_view(view1));
  EXPECT_TRUE(access.isSparse);
  EXPECT_TRUE(access.isRead);
  EXPECT_FALSE(access.isWrite);
  

}

TYPED_TEST_SUITE_P(SparseAccessTest);
template <typename T>
class SparseAccessTest : public ::testing::Test
{
};

TYPED_TEST_P(SparseAccessTest, SymEval)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  SparseAccessTestImpl<EXEC_POLICY>();
 }

REGISTER_TYPED_TEST_SUITE_P(SparseAccessTest, SymEval);
