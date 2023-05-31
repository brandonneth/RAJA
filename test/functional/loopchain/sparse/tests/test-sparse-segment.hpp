#include "RAJA_test-base.hpp"
#include <string>
#include <vector>


template <typename EXEC_POLICY>
void SparseSegmentDenseTestImpl() {
  using namespace RAJA;

  std::stringstream s1;
  for(int i = 0; i < 5; i++) {
    for(int j = 5; j < 10; j++) {
      s1 << i << j << "\n";
    }
  }

  std::stringstream s2;
  auto lam = [&](auto i, auto j) {
    s2 << i << j << "\n";
  };

  using POL=KernelPolicy<
    statement::For<0,loop_exec,
      statement::For<1,loop_exec,
        statement::Lambda<0>
  >>>;
  
  auto seg1 = make_dense_segment(0,5);
  auto seg2 = make_dense_segment(5,10);
  auto segs = make_tuple(seg1,seg2);
  kernel<POL>(segs, lam);
  EXPECT_EQ(s1.str(), s2.str()) << "SparseSegmentDense as segment tuple test\n";
}


template <typename EXEC_POLICY>
void SparseSegmentTestImpl() {
   using namespace RAJA;
   std::vector<camp::idx_t> r{0,1,1,2,2,2,2,3,5,5};
   std::vector<camp::idx_t> c{3,1,4,0,2,3,5,5,2,4};
   std::vector<double>v{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};
  
 
  auto view01 = make_sparse_view<double>(r,c,v);
  std::stringstream s1;
  auto lam = [&](auto i, auto j) {
    s1 << i << ", " << j << " | ";
  };

  std::stringstream correct1;
  for(size_t i = 0; i < r.size(); i++) {
    correct1 << r[i] << ", " << c[i] << " | ";
  }

  using POL=KernelPolicy<
    statement::For<0,loop_exec,
      statement::For<1,loop_exec,
        statement::Lambda<0>
  >>>;
  
  auto seg1_impl = new LeadSegmentImpl(view01.impl.dim(0));
  auto seg1 = SparseSegment(seg1_impl);
  auto seg2 = make_follow_segment(view01.impl.dim(1), seg1_impl);

  auto segs = make_tuple(seg1, seg2);

  kernel<POL>(segs, lam);

  EXPECT_EQ(correct1.str(), s1.str());
   
}

TYPED_TEST_SUITE_P(SparseSegmentTest);
template <typename T>
class SparseSegmentTest : public ::testing::Test
{
};

TYPED_TEST_P(SparseSegmentTest, SparseSegment)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  SparseSegmentDenseTestImpl<EXEC_POLICY>();
  SparseSegmentTestImpl<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SparseSegmentTest, SparseSegment);
