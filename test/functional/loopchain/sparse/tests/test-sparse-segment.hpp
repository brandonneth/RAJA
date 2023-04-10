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

  auto seg0 = SparseSymbolicSegmentLead<2>(view01, 0);
  auto seg1 = SparseSymbolicSegmentFollow<2>(view01, 1, &seg0);

  auto segs = camp::tie(seg0,seg1);

  kernel<POL>(segs, lam);
  EXPECT_EQ(correct1.str(), s1.str());

   std::vector<camp::idx_t> c2{0,1,1,2,2,2,2,3,5,5};
   std::vector<camp::idx_t> r2{3,1,4,0,2,3,5,5,2,4};
   std::vector<double>v2{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};


  auto perm10 = std::vector<camp::idx_t>{1,0};

  std::stringstream s2;

  std::stringstream correct2;
  for(size_t i = 0; i < r2.size(); i++) {
    correct2 << r2[i] << ", " << c2[i] << " | ";
  }


/*
  auto permseg0 = SparseSymbolicSegmentLead<2>(view10, 1);
  auto permseg1 = SparseSymbolicSegmentFollow<2>(view10, 0, &permseg0);

  auto permsegs = camp::tie(permseg0,permseg1);

  s2.str(std::string());
  s2.clear();

  kernel<POL>(permsegs, lam2);

  EXPECT_EQ(correct2.str(), s2.str()) << "Failed for permuted view's sparse iteration space\n";

  auto spaceFromView01 = view01.get_sparse_iteration_space(perm01);

  s1.str(std::string());
  s1.clear();

  kernel<POL>(spaceFromView01, lam);
  EXPECT_EQ(correct1.str(), s1.str()) << "view01's ispace method failed for dimOrder 01\n";
  

  auto spaceFromView10 = view10.get_sparse_iteration_space(perm10);
  
  s1.str(std::string());
  s1.clear();

  kernel<POL>(spaceFromView10, lam);
  EXPECT_EQ(correct1.str(), s1.str()) << "view10's ispace method failed for dimOrder 10\n";
*/
  
}

TYPED_TEST_SUITE_P(SparseSegmentTest);
template <typename T>
class SparseSegmentTest : public ::testing::Test
{
};

TYPED_TEST_P(SparseSegmentTest, SparseSegment)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  SparseSegmentTestImpl<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SparseSegmentTest, SparseSegment);
