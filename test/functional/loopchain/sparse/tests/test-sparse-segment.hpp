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
void CompressedSparseSegmentUnitTestImpl() {
  using namespace RAJA;
  std::vector<idx_t> r{0,1,1,2,2,2,2,3,5,5};
  std::vector<idx_t> c{3,1,4,0,2,3,5,5,2,4};
  std::vector<double>v{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};

  auto compressedLeadSeg = make_compressed_lead_segment(r);
  auto compressedLeadSegImpl = (CompressedLeadSegmentImpl*) compressedLeadSeg.impl;
  
  std::vector<idx_t> correctCompressedRow{0,1,3,7,8,8,10};
  std::vector<idx_t> indexVector = compressedLeadSegImpl->indices;

  auto len = indexVector.size();
  EXPECT_EQ(len, 7) << "compressed lead segment length is not the right length.";
  for(idx_t i = 0; i < correctCompressedRow.size(); i++) {
    auto correctVal = correctCompressedRow[i];
    auto givenVal = indexVector[i];
    EXPECT_EQ(givenVal,correctVal) << "Incorrect index vector in dimension " << i << "\n";
  }

  auto compressedFollowSeg = make_compressed_follow_segment(c, compressedLeadSegImpl);
  
  std::vector<idx_t> rev_c(c.rbegin(), c.rend());

  //simulating the kernel execution of the two dimensional loop.
/*
  //calculating and checking the length of the outer, compressed loop 
  auto spanTuple = make_tuple(compressedLeadSeg.begin(), compressedLeadSeg.end());
  auto innerSpanTuple = make_tuple(compressedFollowSeg.begin(), compressedFollowSeg.end());
  auto leadLen = get<1>(spanTuple) - get<0>(spanTuple);
  EXPECT_EQ(leadLen, correctCompressedRow.size() - 1); //this one is shorter bc its the number of iterations to run, not the number of elements in the row pointer array, which is one longer than the number of iterations the outer loop runs.

  for(int i = 0; i < leadLen; i++) {
    auto leadIteratorValue = get<0>(spanTuple)[i];
    EXPECT_EQ(leadIteratorValue, i);
    
    auto innerStartValue = compressedLeadSegImpl->indices[compressedLeadSegImpl->currIdx];
    auto innerStopValue = compressedLeadSegImpl->indices[compressedLeadSegImpl->currIdx + 1];
    EXPECT_EQ(innerStartValue, correctCompressedRow[i]); 
    EXPECT_EQ(innerStopValue, correctCompressedRow[i+1]); 
   
    EXPECT_EQ(compressedLeadSegImpl->currIdx, i) << "curr index in simulated executoin is not correct\n"; 
    //creating and checking the inner loop values
    auto innerLen = get<1>(innerSpanTuple) - get<0>(innerSpanTuple);
    std::vector<idx_t> correctInnerLens{1, 2, 4, 1, 0, 2};
    EXPECT_EQ(innerLen, correctInnerLens[i]);
    
    for(int j = 0; j < innerLen; j++) {
      auto innerIteratorValue = get<0>(innerSpanTuple)[j];
      
      EXPECT_EQ(innerIteratorValue, rev_c.back()) << "Incorrect inner value for simulated execution with test iterators " << i << " and " << j << "\n";
      rev_c.pop_back();
    }
  }

  compressedLeadSegImpl->currIdx = 0; 
*/





  auto view01 = make_sparse_view<double>(r,c,v);
  std::stringstream s1;
  auto lam = [&](auto i, auto j) {
    //EXPECT_EQ(compressedLeadSegImpl->currIdx, i) << "curr index in lamda execution is not correect\n"; 
    
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
  

  auto segs = make_tuple(compressedLeadSeg, compressedFollowSeg);

  kernel<POL>(segs, lam);

  EXPECT_EQ(correct1.str(), s1.str());

  
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
  CompressedSparseSegmentUnitTestImpl<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SparseSegmentTest, SparseSegment);
