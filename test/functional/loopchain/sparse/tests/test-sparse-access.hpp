#include "RAJA_test-base.hpp"
#include <string>
#include <vector>


template <typename EXEC_POLICY>
void SparseAccessTestSpMV() {
  using namespace RAJA;

  using View1 = RAJA::View<double, RAJA::Layout<1, idx_t, 1>>;

  camp::idx_t Ni = 20;
  camp::idx_t Nj = 30;
  
  SparseView<double, 2> A = make_random_sparse_view2<double>(Ni,Nj, 0,5);
  View1 x(new double[Nj], Nj);
  View1 y(new double[Ni], Ni);

  using POLICY = KernelPolicy<
    statement::For<0, loop_exec,
      statement::For<1, loop_exec,
        statement::Lambda<0>
      >
    >
  >;

  auto seg1 = RangeSegment(0,Ni);
  auto seg2 = RangeSegment(0,Nj);
  auto segs = camp::make_tuple(seg1, seg2);

  auto lam = [&](auto i, auto j) {
    y(i) += A(i,j) * x(j);
  };

  auto knl = make_kernel<POLICY>(segs, lam);
  auto accesses = knl.execute_symbolically();

  for(auto access : accesses) {
    if (access.is_view(x)) {

    }
    if (access.is_view(y)) {

    }
    if (access.is_view(A)) {

    }
  }
}




template <typename EXEC_POLICY>
void SparseAccessTestGauSei() {
  using namespace RAJA;

  using View1 = RAJA::View<double, RAJA::Layout<1, idx_t, 1>>;

  camp::idx_t N = 20;
  
  SparseView<double, 2> A = make_random_sparse_view2<double>(N,N, 0,5);
  View1 b(new double[N], N);
  View1 phi(new double[N], N);
  View1 temp(new double[1], 1);

  using POLICY = KernelPolicy<
    statement::For<0,loop_exec,
      statement::Lambda<0, Segs<0>>,
      statement::For<1,loop_exec,
        statement::Lambda<1>
      >,
      statement::Lambda<2, Segs<0>>
    >
  >;

  auto lam1 = [&](auto) {
    temp(0) = 0.0;
  };
  auto lam2 = [&](auto i, auto j) {
    temp(0) += A(i,j) * phi(j);
  };
  auto lam3 = [&](auto i) {
    phi(i) = (b(i) - temp(0)) / A(i,i);
  };

  auto seg = RangeSegment(0,N);
  auto segs = make_tuple(seg,seg);

  auto knl = make_kernel<POLICY>(segs, lam1, lam2, lam3);
  auto accesses = knl.execute_symbolically();
}

template <typename EXEC_POLICY>
void SparseAccessTestInCholFact() {
  using namespace RAJA;

  camp::idx_t N = 20;
  
  SparseView<double, 2> A = make_random_sparse_view2<double>(N,N, 0,5);
  
  using POLICY = KernelPolicy<
    statement::For<0,loop_exec,
      statement::Lambda<0, Segs<0>>,
      statement::For<1,loop_exec,
        statement::Lambda<1, Segs<0,1>>
      >,
      statement::For<1,loop_exec,
        statement::For<2,loop_exec,
          statement::Lambda<2, Segs<0,1,2>>
        >
      >
    >
  >;

  auto lam1 = [&](auto i0) {
    A(i0,i0) = sqrt(A(i0,i0));
  };
  auto lam2 = [&](auto i0, auto i1) {
    if(A(i1,i0) != 0) {
      A(i1,i0) /= A(i0,i0);
    }
  };
  auto lam3 = [&](auto i0, auto i1, auto i2) {
    if(A(i2,i1) != 0) {
      A(i2,i1) -= A(i2,i0) * A(i1,i0);
    }
  };
  
  auto seg1 = make_symbolic_segment(0,N);
  auto seg2 = make_symbolic_segment(seg1+1, N);
  auto seg3 = make_symbolic_segment(seg2,N);
  auto segs = make_tuple(seg1 ,seg2, seg3);

  auto knl = make_kernel<POLICY>(segs, lam1, lam2, lam3);
  
  auto accesses = knl.execute_symbolically();

}

template <typename EXEC_POLICY>
void SparseAccessTestImpl() {
   using namespace RAJA;
  std::vector<camp::idx_t> r{0,1,1,2,2,2,2,3,5,5};
  std::vector<camp::idx_t> c{3,1,4,0,2,3,5,5,2,4};
  std::vector<double>v{1.0,4.0,9.0,2.0,2.0,3.0,6.0,5.0,7.0,8.0};


  //SparseView<2,double> view1 = make_sparse_view<double>(r,c,v);
  auto view1 = make_sparse_view<double>(r,c,v);

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
  SparseAccessTestSpMV<EXEC_POLICY>();
  SparseAccessTestGauSei<EXEC_POLICY>();
  SparseAccessTestInCholFact<EXEC_POLICY>();
 }

REGISTER_TYPED_TEST_SUITE_P(SparseAccessTest, SymEval);
