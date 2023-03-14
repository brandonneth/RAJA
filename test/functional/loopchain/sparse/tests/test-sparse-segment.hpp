#include "RAJA_test-base.hpp"
#include <string>
#include <vector>
template <typename EXEC_POLICY>
void SparseSegmentTestImpl() {
   using namespace RAJA;
  std::vector<camp::idx_t> r{0,1,1,2,2,2,2,3,5,5};
  std::vector<camp::idx_t> c{3,1,4,0,2,3,5,5,2,4};
  std::vector<double>v{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};


  auto view01 = make_sparse_view<double>(r,c,v);
  auto perm10 = std::vector<camp::idx_t>{1,0};
  //auto view10 = make_sparse_view_permuted<double>(r,c,v,perm10);

  std::stringstream s1;
  auto lam = [&](auto i, auto j) {
    s1 << i << ", " << j << "|";
  };

  std::stringstream correct1;
  for(size_t i = 0; i < r.size(); i++) {
    correct1 << r[i] << ", " << c[i] << "|";
  }

  using POL=KernelPolicy<
    statement::For<0,loop_exec,
      statement::For<1,loop_exec,
        statement::Lambda<0>
  >>>;

  auto seg0 = SparseSymbolicSegmentLead<2>(view01, 0);
  auto seg1 = SparseSymbolicSegmentFollow<2>(view01, 1, &seg0);

  auto segs = camp::tie(seg0,seg1);

  kernel<POL>(segs, lam);
  EXPECT_EQ(correct1.str(), s1.str());
}

TYPED_TEST_SUITE_P(SparseSegmentTest);
template <typename T>
class SparseSegmentTest : public ::testing::Test
{
};

TYPED_TEST_P(SparseSegmentTest, SymEval)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  SparseSegmentTestImpl<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SparseSegmentTest, SymEval);
