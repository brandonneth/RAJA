#include "RAJA_test-base.hpp"
#include <string>
#include <vector>
#include "RAJA/loopchain/Sparse/SparseKernel.hpp"

template <typename EXEC_POLICY>
void CompressedSparseSegmentUnitTestImpl() {
  using namespace RAJA;
  std::vector<idx_t> r{0,1,1,2,2,2,2,3,5,5};
  std::vector<idx_t> c{3,1,4,0,2,3,5,5,2,4};
  std::vector<double>v{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};

  auto compressedFollowSeg = make_compressed_follow_segment(c);
  

  auto compressedLeadSeg = make_compressed_lead_segment(r, compressedFollowSeg);
  
  std::vector<idx_t> correctCompressedRow{0,1,3,7,8,8,10};
  
  auto compressedBegin = (CompressedLeadIterator*)compressedLeadSeg.begin().impl;

  for(size_t i = 0; i < correctCompressedRow.size(); i++) {
    EXPECT_EQ(compressedBegin->dimStarts[i], correctCompressedRow[i]);
  }
  
  auto view01 = make_sparse_view<double>(r,c,v);
  std::stringstream s1;
  auto lam = [&](auto i, auto j) {
    double temp = view01(i,j); 
    camp::sink(temp);
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

  auto dense_segs = make_tuple(RangeSegment(0,6), RangeSegment(0,6));

  auto dense_knl = make_kernel<POL>(dense_segs, lam);
  auto accesses = dense_knl.template execute_symbolically();
  EXPECT_EQ(accesses.size(), 1);

  auto sknl = make_sparse_kernel_sym<POL>(make_tuple(0,1), dense_segs, view01, lam);

  s1.str("");
  sknl();
  
  EXPECT_EQ(correct1.str(), s1.str()) << "perfect nesting fails sparsification when using make_sparse_kernel for compressed dims\n";
 

}

template <typename EXEC_POLICY>
void SparseSegmentGauSei() {
  using namespace RAJA;

  size_t N = 5;
  auto view = make_random_sparse_view2<double,true>(N,1.0);
  std::stringstream s2;

  using POLICY = KernelPolicy<
      statement::For<0,loop_exec,
        statement::Lambda<0,Segs<0>>,
        statement::For<1,loop_exec,
          statement::Lambda<1>
        >,
        statement::Lambda<2,Segs<0>>
      >
    >;

  auto lam1 = [&](auto i) {
    s2 << "1" << i << " ";
  };
  auto lam2 = [&](auto i, auto j) {
    s2 << "2" << i << j << " ";
    if (j != i) {
      auto temp = 0.0;
      temp += view(i,j);
    }
  };
  auto lam3 = [&](auto i) {
    s2 << "3" << i << " "; 
  };

  auto seg1 = RangeSegment(0,N);
  auto seg2 = RangeSegment(0,N);
  auto dense_segs = make_tuple(seg1, seg2);
 
  auto dense_gau = make_kernel<POLICY>(dense_segs, lam1, lam2, lam3);
  s2.str("");
  dense_gau();
  auto correct2 = s2.str();
   
  for(size_t i = 0; i < N*N; i++) {
    EXPECT_EQ(view.dim(0)[i], i / 5) << "pre failed for index " << i << "\n";
    EXPECT_EQ(view.dim(1)[i], i % 5)<< "pre failed for index " << i << "\n";

  }

  auto sparse_gau = make_sparse_kernel_sym<POLICY,1>(make_tuple(0,1), dense_segs, view, lam1, lam2, lam3);
  
  for(size_t i = 0; i < N*N; i++) {
    EXPECT_EQ(view.dim(0)[i], i / 5) << "failed for index " << i << "\n";
    EXPECT_EQ(view.dim(1)[i], i % 5)<< "failed for index " << i << "\n";

  }

  s2.str("");
  sparse_gau();
  auto output2 = s2.str();

  EXPECT_STREQ(output2.c_str(), correct2.c_str()) << "Failed schedule for gausei style loop.\n";
}


template <typename EXEC_POLICY>
void CompressedSegment3DTest() {
  using namespace RAJA;

  size_t N = 5;
  auto view = make_random_sparse_view3<double>(N);
  std::stringstream s2;
  double temp; 
  auto lam = [&](auto i, auto j, auto k) {
    s2 << "lam" << i << j << k << " ";
    double temp;
    temp = view(i,j,k);
  };

  using POLICY=KernelPolicy<
    statement::For<0,loop_exec,
    statement::For<1,loop_exec,
    statement::For<2,loop_exec,
    statement::Lambda<0>
  >>>>;

  std::stringstream correct;
  auto correct_lam = [&](auto i, auto j, auto k) {
    if (view(i,j,k) != 0) {
      correct << "lam" << i << j << k << " ";
    }
  };
  auto seg1 = RangeSegment(0,N);
  auto dense_segs = make_tuple(seg1, seg1, seg1);
  kernel<POLICY>(dense_segs, correct_lam);

  auto sparse_knl = make_sparse_kernel_sym<POLICY>(make_tuple(0,0,0), dense_segs, view, lam);
  s2.str(""); 
  sparse_knl();

  EXPECT_EQ(s2.str(), correct.str());

}




template <typename ArgTuple, typename Lambda, camp::idx_t...Is>
auto wrap_sparse_lambda(const Lambda & l,camp::idx_seq<Is...>) {
  return [=](ArgTuple tuple) {
    l(camp::get<Is>(tuple)...);
  };
}




template <typename EXEC_POLICY>
void IndexSetTest() {
  using namespace RAJA;
  std::stringstream out;

  auto twoArgs = [&](auto i, auto j) {
    out << "twoArgs" << " " << i << " " << j << "\n";
  };

  auto tupleArg = [&](auto tup) {
    auto i = camp::get<0>(tup);
    auto j = camp::get<1>(tup);
    out << "tuple" << " " << i << " " << j << "\n";
  };

  auto wrapped = wrap_sparse_lambda<camp::tuple<int,int>>(twoArgs,camp::idx_seq<0,1>{});

  using TupleVector = std::vector<camp::tuple<int,int>>;
  TupleVector tuples; 
  tuples.push_back(make_tuple(0,0));
  tuples.push_back(make_tuple(1,1));
  tuples.push_back(make_tuple(1,8));
  tuples.push_back(make_tuple(3,2));
  tuples.push_back(make_tuple(3,3));

  TypedIndexSet<TupleVector> iset;
  iset.push_back(tuples);

  using IndexPolicy = ExecPolicy<seq_segit, seq_exec>;
  
  forall<IndexPolicy>(iset, tupleArg);

  std::string expected("tuple 0 0\ntuple 1 1\ntuple 1 8\ntuple 3 2\ntuple 3 3\n");

  EXPECT_EQ(out.str(), expected);

  out.str("");
  out.clear();

  forall<IndexPolicy>(iset, wrapped);

  std::string expected2("twoArgs 0 0\ntwoArgs 1 1\ntwoArgs 1 8\ntwoArgs 3 2\ntwoArgs 3 3\n");
  
  EXPECT_EQ(out.str(), expected2) << "wrapped index tuple lambda didn't work.\n";
};




TYPED_TEST_SUITE_P(SparseSegmentTest);
template <typename T>
class SparseSegmentTest : public ::testing::Test
{
};



TYPED_TEST_P(SparseSegmentTest, SparseSegment)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  CompressedSparseSegmentUnitTestImpl<EXEC_POLICY>();
  SparseSegmentGauSei<EXEC_POLICY>() ;
  SparseSegmentGauSei<EXEC_POLICY>() ;
  IndexSetTest<EXEC_POLICY>();
  CompressedSegment3DTest<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(SparseSegmentTest, SparseSegment);
