#include "RAJA_test-base.hpp"
#include <string>
#include <vector>
#include <bitset>
#include "RAJA/loopchain/Sparse/SparseView.hpp"





template <typename EXEC_POLICY>
void SparseViewTestImpl() {
   using namespace RAJA;
  std::vector<camp::idx_t> r{0,1,1,2,2,2,2,3,5,5};
  std::vector<camp::idx_t> c{3,1,4,0,2,3,5,5,2,4};
  std::vector<double>v{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};


  //SparseView<2,double> view1 = make_sparse_view<double>(r,c,v);
  auto view1 = make_sparse_view<double,2,COO,std::vector<idx_t>,std::vector<idx_t>,std::vector<double>>(r,c,v);
  for(size_t i = 0; i < r.size(); i++ ) {
    EXPECT_EQ(view1(r[i], c[i]), v[i]) << "Failed for i = " << i << ".";
  }

  auto perm10 = std::vector<camp::idx_t>{1,0};
  auto view10 = make_sparse_view_permuted<double,2,COO>(r,c,v,perm10);

  //making sure the original vectors aren't being touched.
  EXPECT_EQ(r[0], 0);
  EXPECT_EQ(r[2], 1);
  EXPECT_EQ(r[4], 2);
  EXPECT_EQ(r[6], 2);
  EXPECT_EQ(r[8], 5);
  EXPECT_EQ(c[0], 3);
  EXPECT_EQ(c[2], 4);
  EXPECT_EQ(c[4], 2);
  EXPECT_EQ(c[6], 5);
  EXPECT_EQ(c[8], 2);
  

  //make sure the permuted view's indices are sorted by the column values.
  EXPECT_EQ(view10.dim(1)[0], 0);
  EXPECT_EQ(view10.dim(1)[1], 1);
  EXPECT_EQ(view10.dim(1)[2], 2);
  EXPECT_EQ(view10.dim(1)[3], 2);
  EXPECT_EQ(view10.dim(1)[4], 3);
  EXPECT_EQ(view10.dim(1)[5], 3);
  EXPECT_EQ(view10.dim(1)[6], 4);
  EXPECT_EQ(view10.dim(1)[7], 4);
  EXPECT_EQ(view10.dim(1)[8], 5);
  EXPECT_EQ(view10.dim(1)[9], 5);
  
  EXPECT_EQ(view10.dim(0)[0], 2);
  EXPECT_EQ(view10.dim(0)[1], 1);
  EXPECT_EQ(view10.dim(0)[2], 2);
  EXPECT_EQ(view10.dim(0)[3], 5);
  EXPECT_EQ(view10.dim(0)[4], 0);
  EXPECT_EQ(view10.dim(0)[5], 2);
  EXPECT_EQ(view10.dim(0)[6], 1);
  EXPECT_EQ(view10.dim(0)[7], 5);
  EXPECT_EQ(view10.dim(0)[8], 2);
  EXPECT_EQ(view10.dim(0)[9], 3);

  EXPECT_EQ(view10.impl.permutation[0], 1);
  EXPECT_EQ(view10.impl.permutation[1], 0);
  EXPECT_EQ(view10.val(0), 2.0) << "First value in permuted sparse view not correct\n";
  EXPECT_EQ(view10.val(1), 4.0) << "Second value in permuted sparse view not correct\n";
  EXPECT_EQ(view10.val(9), 5.0) << "Last value in permuted sparse view not correct\n";
  EXPECT_EQ(1.0, view10(0,3)) << "Permuted sparse view fails access on (0,3)";
  EXPECT_EQ(2.0, view10(2,0)) << "Permuted sparse view fails access on (2,0)";
  EXPECT_EQ(0.0, view10(3,0)) << "Permuted sparse view fails access on (3,0)";

  std::array<idx_t, 2> findArg{{0,3}};
  auto findIdx = view10.impl.find<0>(0,10, findArg);
  EXPECT_EQ(findIdx,4); 

}

template <typename EXEC_POLICY>
void SparseViewDiagTest() {
  using namespace RAJA;
  std::vector<camp::idx_t> r{0,1,1,2,2,2,2,3,4,5};
  std::vector<camp::idx_t> c{0,1,2,0,2,3,5,3,4,5};
  std::vector<double>v{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};

  auto view = make_sparse_view<double,2,DIAG>(r,c,v);
  
  auto d = view(0,0);
  d = view(1,1);
  EXPECT_EQ(view.get_hit_rate(), 1.0) << "DIAG view not hitting on diagonal entries 00 or 11";
  d = view(2,2);
  EXPECT_EQ(view.get_hit_rate(), 1.0) << "DIAG view not hitting on diagonal entries 22";
  d = view(3,3);
  EXPECT_EQ(view.get_hit_rate(), 1.0) << "DIAG view not hitting on diagonal entries 22";

}

template <typename EXEC_POLICY>
void SparseViewCSRTest() {
  using namespace RAJA;
  std::vector<camp::idx_t> r{0,1,1,2,2,2,2,3,4,5};
  std::vector<camp::idx_t> c{0,1,2,0,2,3,5,3,4,5};
  std::vector<double>v{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};

  auto view = make_sparse_view<double,2,CSR>(r,c,v);
  
  auto d = view(0,0);
  EXPECT_EQ(view.get_hit_rate(), 1.0) << "CSR fails expected access on first value.";
  d = view(1,1);
  EXPECT_EQ(view.get_hit_rate(), 1.0) << "CSR fails expected access on 1 1.";
  d = view(1,2);
  EXPECT_EQ(view.get_hit_rate(), 1.0) << "CSR fails expected access on 1 2.";
  d = view(2,0);
  EXPECT_EQ(view.get_hit_rate(), 1.0) << "CSR fails expected access on 2 0.";
  d = view(1,1);
  EXPECT_EQ(view.get_hit_rate(), 0.8) << "CSR does not fail expected access as expectde.";
 
  std::vector<idx_t> perm10{{1,0}};
  auto view2 = make_sparse_view_permuted<double,2,CSR>(r,c,v,perm10);
   d = view2(0,0);
  EXPECT_EQ(view2.get_hit_rate(), 1.0) << "CSC fails expected access on first value.";
  d = view2(2,0);
  EXPECT_EQ(view2.get_hit_rate(), 1.0) << "CSC fails expected access on 2 0.";
  d = view2(1,1);
  EXPECT_EQ(view2.get_hit_rate(), 1.0) << "CSC fails expected access on 1,1.";
  d = view2(1,2);
  EXPECT_EQ(view2.get_hit_rate(), 1.0) << "CSC fails expected access on 2 0.";
  d = view2(5,5);
  EXPECT_EQ(view2.get_hit_rate(), 0.8) << "CSC does not fail expected access as expectde." << d;

}


template <typename EXEC_POLICY>
void make_sparse_view_test() {
  using namespace RAJA;
  std::vector<camp::idx_t> r{0,1,1,2,2,2,2,3,4,5};
  std::vector<camp::idx_t> c{0,1,2,0,2,3,5,3,4,5};
  std::vector<double>v{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};

  auto viewCOO = make_sparse_view<double, 2, COO>(r,c,v);
  auto viewDIAG = make_sparse_view<double,2, DIAG>(r,c,v);
  auto viewCSR = make_sparse_view<double,2, CSR>(r,c,v);

  std::vector<idx_t> perm10{{1,0}};
  auto viewCOO_10 = make_sparse_view_permuted<double, 2, COO>(r,c,v,perm10);
  auto viewDIAG_10 = make_sparse_view_permuted<double, 2, DIAG>(r,c,v,perm10);

  std::string cooname = COO<double, 2>::name;
  std::string csrname = CSR<double, 2>::name;
  std::string diagname = DIAG<double, 2>::name;
  EXPECT_EQ(viewDIAG_10.impl.fill, viewDIAG.impl.fill);
  EXPECT_EQ(viewDIAG_10.impl.fill, viewCOO_10.impl.fill);
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      EXPECT_EQ(viewCOO(i,j),viewDIAG(i,j));
      auto c10 = viewCOO_10(i,j);
      auto d10 = viewDIAG_10(i,j);
      if (c10 != d10) {
        std::cout << "failed conditional check\n";
      }
      EXPECT_EQ(c10,d10) << "iteration: " << i << " " << j << "\n " 
      << "coo bits: " << std::bitset<64>(c10) << "\n"
      << "dia bits: " << std::bitset<64>(d10) << "\n";
     EXPECT_EQ(c10,d10);
      EXPECT_EQ(viewCOO_10(i,j),viewDIAG_10(i,j)) << "iteration: " << i << " " << j << "\n " 
      << "coo bits: " << std::bitset<8>(viewCOO_10(i,j)) << "\n"
      << "dia bits: " << std::bitset<8>(viewDIAG_10(i,j)) << "\n";
      EXPECT_EQ(viewCOO(i,j),viewDIAG_10(i,j));
      EXPECT_EQ(viewCOO(i,j),viewCSR(i,j));
    }
  }
  for(size_t i = 0; i < r.size(); i++) {
    EXPECT_EQ(viewCOO(r[i],c[i]), v[i]);
    EXPECT_EQ(viewCOO_10(r[i],c[i]), v[i]);
    EXPECT_EQ(viewDIAG(r[i],c[i]), v[i]);
    EXPECT_EQ(viewDIAG_10(r[i],c[i]), v[i]);
    EXPECT_EQ(viewCSR(r[i],c[i]), v[i]);
  }  
}



TYPED_TEST_SUITE_P(SparseViewTest);
template <typename T>
class SparseViewTest : public ::testing::Test
{
};

TYPED_TEST_P(SparseViewTest, SymEval)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  SparseViewTestImpl<EXEC_POLICY>();
  SparseViewDiagTest<EXEC_POLICY>();
  SparseViewCSRTest<EXEC_POLICY>();
  make_sparse_view_test<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SparseViewTest, SymEval);
