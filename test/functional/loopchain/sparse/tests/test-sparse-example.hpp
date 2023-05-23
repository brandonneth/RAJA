#include "RAJA_test-base.hpp"
#include <string>
#include <vector>


template <typename T>
void print_perm(T A) {
  std::cout << "PERM: ";
  for(auto i : A.impl->get_permutation()) {
    std::cout << i << " ";
  }
  std::cout << "\n";
}

template <typename EXEC_POLICY>
void SparseExampleTest3Impl() {
  using namespace RAJA;

  camp::idx_t N = 3;

  std::vector<idx_t> d0 = {0,0,1,1,2,2,0};
  std::vector<idx_t> d1 = {1,2,1,2,0,0,2};
  std::vector<idx_t> d2 = {0,0,0,1,1,2,2};
  std::vector<double> v = {1,2,3,4,5,6,7};
  auto A = make_sparse_view<double>(d0,d1,d2,v);

  using POL012 =  KernelPolicy<statement::For<0,loop_exec,statement::For<1,loop_exec,statement::For<2,loop_exec,statement::Lambda<0>>>>>;
  using POL021 =  KernelPolicy<statement::For<0,loop_exec,statement::For<2,loop_exec,statement::For<1,loop_exec,statement::Lambda<0>>>>>;
  using POL201 =  KernelPolicy<statement::For<2,loop_exec,statement::For<0,loop_exec,statement::For<1,loop_exec,statement::Lambda<0>>>>>;

  auto lam012 = [&](auto i0, auto i1, auto i2) {
    if (A(i0,i1,i2) != 0.0) {std::cout << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam021 = [&](auto i0, auto i1, auto i2) {
    if (A(i0,i2,i1) != 0.0) {std::cout << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam102 = [&](auto i0, auto i1, auto i2) {
    if (A(i1,i0,i2) != 0.0) {std::cout << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam120 = [&](auto i0, auto i1, auto i2) {
    if (A(i1,i2,i0) != 0.0) {std::cout << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam201 = [&](auto i0, auto i1, auto i2) {
    if (A(i2,i0,i1) != 0.0) {std::cout << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam210 = [&](auto i0, auto i1, auto i2) {
    if (A(i2,i1,i0) != 0.0) {std::cout << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam_prt = [&](auto i0, auto i1, auto i2) {
    std::cout << i0 << "," << i1 << "," << i2 << "  ";
  };
  auto dense_segs = tuple_repeat<3>(RangeSegment(0,N));

  print_perm(A);
  std::cout << "012 dns 012: ";
  kernel<POL012>(dense_segs, lam012);
  std::cout << "\n";

  std::cout << "201 dns 012: ";
  kernel<POL201>(dense_segs, lam012);
  std::cout << "\n";

  std::cout << "012 dns 201: ";
  kernel<POL012>(dense_segs, lam201);
  std::cout << "\n";
 
  std::vector<camp::idx_t> perm201 = {2,0,1};
  A.impl->set_permutation(perm201);
  auto seg201_0 = SparseSymbolicSegmentLead<3>(A,2);
  auto seg201_1 = SparseSymbolicSegmentFollow<3>(A,0, &seg201_0);
  auto seg201_2 = SparseSymbolicSegmentFollow<3>(A,1, &seg201_0);
  auto segs201 = camp::tie(seg201_0, seg201_1,seg201_2);

  print_perm(A);
  std::cout << "012 201 201: ";
  kernel<POL012>(segs201, lam201);
  std::cout << "\n";


  std::vector<camp::idx_t> perm120 = {1,2,0};
  A.impl->set_permutation(perm120);
  auto seg120_0 = SparseSymbolicSegmentLead<3>(A,1);
  auto seg120_1 = SparseSymbolicSegmentFollow<3>(A,2, &seg120_0);
  auto seg120_2 = SparseSymbolicSegmentFollow<3>(A,0, &seg120_0);
  auto segs120 = camp::tie(seg120_0, seg120_1,seg120_2);
 
  print_perm(A);

  std::cout << "012 120 201: ";
  kernel<POL012>(segs120, lam201);
  std::cout << "\n";

  std::cout << "021 dns 201: ";
  kernel<POL021>(dense_segs, lam201);
  std::cout << "\n";

  std::cout << "021 120 201: ";
  kernel<POL021>(dense_segs, lam201);
  std::cout << "\n";
}


template <typename EXEC_POLICY>
void SparseExampleTestImpl() {
  using namespace RAJA;
  using View1 = RAJA::View<double, RAJA::Layout<1, idx_t, 1>>;

  camp::idx_t Ni = 10;
  camp::idx_t Nj = 10;
  
  SparseView<double, 2> A = make_random_sparse_view2<double>(Ni,Nj, 0,5);
  
  View1 x(new double[Nj], Nj);
  View1 y(new double[Ni], Ni);
  View1 z(new double[Nj], Nj);
  for(int i = 0; i < Nj; i++) {
    x(i) = i * 0.75;
    z(i) = 0.0;
  }
  for(int i = 0; i < Ni; i++) {
    y(i) = 0.0;
  }
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

  auto lam_ij = [&](auto i, auto j) {
    if(A(i,j) != 0) {
      std::cout << i << "," << j << "  ";
    }
  };

  auto knl1 = make_kernel<POLICY>(segs, lam_ij);


  auto lam_ji = [&](auto i, auto j) {
     if(A(j,i) != 0) {
      std::cout << i << "," << j << "  ";
    }  };

  auto segs2 = camp::make_tuple(seg1, seg2);
  auto knl2 = make_kernel<POLICY>(segs2, lam_ji);

  knl1();
  std::cout << "\n";
  knl2();

  auto permseg0 = SparseSymbolicSegmentLead<2>(A, 0);
  auto permseg1 = SparseSymbolicSegmentFollow<2>(A, 1, &permseg0);
  
  auto permsegs = camp::tie(permseg0, permseg1);

  std::cout << "\n";
  kernel<POLICY>(permsegs, [=](auto i, auto j) {
    std::cout << i << "," << j << "  ";
  });

  std::vector<camp::idx_t> perm10 = {1,0};
  A.impl->set_permutation(perm10);

  auto pseg0 = SparseSymbolicSegmentLead<2>(A,1);
  auto pseg1 = SparseSymbolicSegmentFollow<2>(A,0, &pseg0);
  auto psegs = camp::tie(pseg0, pseg1);

  std::cout << "\n";
  kernel<POLICY>(psegs, lam_ji);
  std::cout << "\n";
}

template <typename EXEC_POLICY>
void SparseExampleTestSparsify3() {
  using namespace RAJA;
  
  camp::idx_t N = 10;
  auto A = make_random_sparse_view3<double>(N);
  auto A2 = make_random_sparse_view2<double>(N);
  auto A4 = make_random_sparse_view4<double>(N);
  std::stringstream stream;

  auto lam012 = [&](auto i0, auto i1, auto i2) {
    if (A(i0,i1,i2) != 0.0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam021 = [&](auto i0, auto i1, auto i2) {
    if (A(i0,i2,i1) != 0.0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam102 = [&](auto i0, auto i1, auto i2) {
    if (A(i1,i0,i2) != 0.0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam120 = [&](auto i0, auto i1, auto i2) {
    if (A(i1,i2,i0) != 0.0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam201 = [&](auto i0, auto i1, auto i2) {
    if (A(i2,i0,i1) != 0.0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam210 = [&](auto i0, auto i1, auto i2) {
    if (A(i2,i1,i0) != 0.0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam_prt = [&](auto i0, auto i1, auto i2) {
    stream << i0 << "," << i1 << "," << i2 << "  ";
  };

  auto dense_segs = tuple_repeat<3>(RangeSegment(0,N));
  
  std::vector<camp::idx_t> perm012 = {0,1,2};
  std::vector<camp::idx_t> perm021 = {0,2,1};
  std::vector<camp::idx_t> perm102 = {1,0,2};
  std::vector<camp::idx_t> perm120 = {1,2,0};
  std::vector<camp::idx_t> perm201 = {2,0,1};
  std::vector<camp::idx_t> perm210 = {2,1,0};
 
  using POL012 =  KernelPolicy<statement::For<0,loop_exec,statement::For<1,loop_exec,statement::For<2,loop_exec,statement::Lambda<0>>>>>;
  using POL021 =  KernelPolicy<statement::For<0,loop_exec,statement::For<2,loop_exec,statement::For<1,loop_exec,statement::Lambda<0>>>>>;
  using POL102 =  KernelPolicy<statement::For<1,loop_exec,statement::For<0,loop_exec,statement::For<2,loop_exec,statement::Lambda<0>>>>>;
  using POL120 =  KernelPolicy<statement::For<1,loop_exec,statement::For<2,loop_exec,statement::For<0,loop_exec,statement::Lambda<0>>>>>;
  using POL201 =  KernelPolicy<statement::For<2,loop_exec,statement::For<0,loop_exec,statement::For<1,loop_exec,statement::Lambda<0>>>>>;
  using POL210 =  KernelPolicy<statement::For<2,loop_exec,statement::For<1,loop_exec,statement::For<0,loop_exec,statement::Lambda<0>>>>>;

  
  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam012);
  std::string ref_output_012_012 = stream.str();
  auto sparsified_012_012 = sparsify_iteration_space<POL012>(dense_segs, A, perm012);
  stream.clear();
  stream.str("");
  kernel<POL012>(sparsified_012_012, lam_prt);
  std::string test_output_012_012 = stream.str();
  EXPECT_EQ(test_output_012_012,ref_output_012_012) << "Failed sparsification with 012->012 nesting order change and 012->012 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";
  
  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam021);
  std::string ref_output_012_021 = stream.str();
  auto sparsified_012_021 = sparsify_iteration_space<POL012>(dense_segs, A, perm021);
  stream.clear();
  stream.str("");
  kernel<POL012>(sparsified_012_021, lam_prt);
  std::string test_output_012_021 = stream.str();
  EXPECT_EQ(test_output_012_021,ref_output_012_021) << "Failed sparsification with 012->012 nesting order change and 012->021 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam102);
  std::string ref_output_012_102 = stream.str();
  auto sparsified_012_102 = sparsify_iteration_space<POL012>(dense_segs, A, perm102);
  stream.clear();
  stream.str("");
  kernel<POL012>(sparsified_012_102, lam_prt);
  std::string test_output_012_102 = stream.str();
  EXPECT_EQ(test_output_012_102,ref_output_012_102) << "Failed sparsification with 012->012 nesting order change and 012->102 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam120);
  std::string ref_output_012_120 = stream.str();
  auto sparsified_012_120 = sparsify_iteration_space<POL012>(dense_segs, A, perm120);
  stream.clear();
  stream.str("");
  kernel<POL012>(sparsified_012_120, lam_prt);
  std::string test_output_012_120 = stream.str();
  EXPECT_EQ(test_output_012_120,ref_output_012_120) << "Failed sparsification with 012->012 nesting order change and 012->120 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam201);
  std::string ref_output_012_201 = stream.str();
  auto sparsified_012_201 = sparsify_iteration_space<POL012>(dense_segs, A, perm201);
  stream.clear();
  stream.str("");
  kernel<POL012>(sparsified_012_201, lam_prt);
  std::string test_output_012_201 = stream.str();
  EXPECT_EQ(test_output_012_201,ref_output_012_201) << "Failed sparsification with 012->012 nesting order change and 012->201 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam210);
  std::string ref_output_012_210 = stream.str();
  auto sparsified_012_210 = sparsify_iteration_space<POL012>(dense_segs, A, perm210);
  stream.clear();
  stream.str("");
  kernel<POL012>(sparsified_012_210, lam_prt);
  std::string test_output_012_210 = stream.str();
  EXPECT_EQ(test_output_012_210,ref_output_012_210) << "Failed sparsification with 012->012 nesting order change and 012->210 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL021>(dense_segs, lam012);
  std::string ref_output_021_012 = stream.str();
  auto sparsified_021_012 = sparsify_iteration_space<POL021>(dense_segs, A, perm012);
  stream.clear();
  stream.str("");
  kernel<POL021>(sparsified_021_012, lam_prt);
  std::string test_output_021_012 = stream.str();
  EXPECT_EQ(test_output_021_012,ref_output_021_012) << "Failed sparsification with 012->021 nesting order change and 012->012 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL021>(dense_segs, lam021);
  std::string ref_output_021_021 = stream.str();
  auto sparsified_021_021 = sparsify_iteration_space<POL021>(dense_segs, A, perm021);
  stream.clear();
  stream.str("");
  kernel<POL021>(sparsified_021_021, lam_prt);
  std::string test_output_021_021 = stream.str();
  EXPECT_EQ(test_output_021_021,ref_output_021_021) << "Failed sparsification with 012->021 nesting order change and 012->021 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL021>(dense_segs, lam102);
  std::string ref_output_021_102 = stream.str();
  auto sparsified_021_102 = sparsify_iteration_space<POL021>(dense_segs, A, perm102);
  stream.clear();
  stream.str("");
  kernel<POL021>(sparsified_021_102, lam_prt);
  std::string test_output_021_102 = stream.str();
  EXPECT_EQ(test_output_021_102,ref_output_021_102) << "Failed sparsification with 012->021 nesting order change and 012->102 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL021>(dense_segs, lam120);
  std::string ref_output_021_120 = stream.str();
  auto sparsified_021_120 = sparsify_iteration_space<POL021>(dense_segs, A, perm120);
  stream.clear();
  stream.str("");
  kernel<POL021>(sparsified_021_120, lam_prt);
  std::string test_output_021_120 = stream.str();
  EXPECT_EQ(test_output_021_120,ref_output_021_120) << "Failed sparsification with 012->021 nesting order change and 012->120 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL021>(dense_segs, lam201);
  std::string ref_output_021_201 = stream.str();
  auto sparsified_021_201 = sparsify_iteration_space<POL021>(dense_segs, A, perm201);
  stream.clear();
  stream.str("");
  kernel<POL021>(sparsified_021_201, lam_prt);
  std::string test_output_021_201 = stream.str();
  EXPECT_EQ(test_output_021_201,ref_output_021_201) << "Failed sparsification with 012->021 nesting order change and 012->201 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL021>(dense_segs, lam210);
  std::string ref_output_021_210 = stream.str();
  auto sparsified_021_210 = sparsify_iteration_space<POL021>(dense_segs, A, perm210);
  stream.clear();
  stream.str("");
  kernel<POL021>(sparsified_021_210, lam_prt);
  std::string test_output_021_210 = stream.str();
  EXPECT_EQ(test_output_021_210,ref_output_021_210) << "Failed sparsification with 012->021 nesting order change and 012->210 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(dense_segs, lam012);
  std::string ref_output_102_012 = stream.str();
  auto sparsified_102_012 = sparsify_iteration_space<POL102>(dense_segs, A, perm012);
  stream.clear();
  stream.str("");
  kernel<POL102>(sparsified_102_012, lam_prt);
  std::string test_output_102_012 = stream.str();
  EXPECT_EQ(test_output_102_012,ref_output_102_012) << "Failed sparsification with 012->102 nesting order change and 012->012 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(dense_segs, lam021);
  std::string ref_output_102_021 = stream.str();
  auto sparsified_102_021 = sparsify_iteration_space<POL102>(dense_segs, A, perm021);
  stream.clear();
  stream.str("");
  kernel<POL102>(sparsified_102_021, lam_prt);
  std::string test_output_102_021 = stream.str();
  EXPECT_EQ(test_output_102_021,ref_output_102_021) << "Failed sparsification with 012->102 nesting order change and 012->021 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(dense_segs, lam102);
  std::string ref_output_102_102 = stream.str();
  auto sparsified_102_102 = sparsify_iteration_space<POL102>(dense_segs, A, perm102);
  stream.clear();
  stream.str("");
  kernel<POL102>(sparsified_102_102, lam_prt);
  std::string test_output_102_102 = stream.str();
  EXPECT_EQ(test_output_102_102,ref_output_102_102) << "Failed sparsification with 012->102 nesting order change and 012->102 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(dense_segs, lam120);
  std::string ref_output_102_120 = stream.str();
  auto sparsified_102_120 = sparsify_iteration_space<POL102>(dense_segs, A, perm120);
  stream.clear();
  stream.str("");
  kernel<POL102>(sparsified_102_120, lam_prt);
  std::string test_output_102_120 = stream.str();
  EXPECT_EQ(test_output_102_120,ref_output_102_120) << "Failed sparsification with 012->102 nesting order change and 012->120 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(dense_segs, lam201);
  std::string ref_output_102_201 = stream.str();
  auto sparsified_102_201 = sparsify_iteration_space<POL102>(dense_segs, A, perm201);
  stream.clear();
  stream.str("");
  kernel<POL102>(sparsified_102_201, lam_prt);
  std::string test_output_102_201 = stream.str();
  EXPECT_EQ(test_output_102_201,ref_output_102_201) << "Failed sparsification with 012->102 nesting order change and 012->201 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(dense_segs, lam210);
  std::string ref_output_102_210 = stream.str();
  auto sparsified_102_210 = sparsify_iteration_space<POL102>(dense_segs, A, perm210);
  stream.clear();
  stream.str("");
  kernel<POL102>(sparsified_102_210, lam_prt);
  std::string test_output_102_210 = stream.str();
  EXPECT_EQ(test_output_102_210,ref_output_102_210) << "Failed sparsification with 012->102 nesting order change and 012->210 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL120>(dense_segs, lam012);
  std::string ref_output_120_012 = stream.str();
  auto sparsified_120_012 = sparsify_iteration_space<POL120>(dense_segs, A, perm012);
  stream.clear();
  stream.str("");
  kernel<POL120>(sparsified_120_012, lam_prt);
  std::string test_output_120_012 = stream.str();
  EXPECT_EQ(test_output_120_012,ref_output_120_012) << "Failed sparsification with 012->120 nesting order change and 012->012 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL120>(dense_segs, lam021);
  std::string ref_output_120_021 = stream.str();
  auto sparsified_120_021 = sparsify_iteration_space<POL120>(dense_segs, A, perm021);
  stream.clear();
  stream.str("");
  kernel<POL120>(sparsified_120_021, lam_prt);
  std::string test_output_120_021 = stream.str();
  EXPECT_EQ(test_output_120_021,ref_output_120_021) << "Failed sparsification with 012->120 nesting order change and 012->021 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL120>(dense_segs, lam102);
  std::string ref_output_120_102 = stream.str();
  auto sparsified_120_102 = sparsify_iteration_space<POL120>(dense_segs, A, perm102);
  stream.clear();
  stream.str("");
  kernel<POL120>(sparsified_120_102, lam_prt);
  std::string test_output_120_102 = stream.str();
  EXPECT_EQ(test_output_120_102,ref_output_120_102) << "Failed sparsification with 012->120 nesting order change and 012->102 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL120>(dense_segs, lam120);
  std::string ref_output_120_120 = stream.str();
  auto sparsified_120_120 = sparsify_iteration_space<POL120>(dense_segs, A, perm120);
  stream.clear();
  stream.str("");
  kernel<POL120>(sparsified_120_120, lam_prt);
  std::string test_output_120_120 = stream.str();
  EXPECT_EQ(test_output_120_120,ref_output_120_120) << "Failed sparsification with 012->120 nesting order change and 012->120 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL120>(dense_segs, lam201);
  std::string ref_output_120_201 = stream.str();
  auto sparsified_120_201 = sparsify_iteration_space<POL120>(dense_segs, A, perm201);
  stream.clear();
  stream.str("");
  kernel<POL120>(sparsified_120_201, lam_prt);
  std::string test_output_120_201 = stream.str();
  EXPECT_EQ(test_output_120_201,ref_output_120_201) << "Failed sparsification with 012->120 nesting order change and 012->201 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL120>(dense_segs, lam210);
  std::string ref_output_120_210 = stream.str();
  auto sparsified_120_210 = sparsify_iteration_space<POL120>(dense_segs, A, perm210);
  stream.clear();
  stream.str("");
  kernel<POL120>(sparsified_120_210, lam_prt);
  std::string test_output_120_210 = stream.str();
  EXPECT_EQ(test_output_120_210,ref_output_120_210) << "Failed sparsification with 012->120 nesting order change and 012->210 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(dense_segs, lam012);
  std::string ref_output_201_012 = stream.str();
  auto sparsified_201_012 = sparsify_iteration_space<POL201>(dense_segs, A, perm012);
  stream.clear();
  stream.str("");
  kernel<POL201>(sparsified_201_012, lam_prt);
  std::string test_output_201_012 = stream.str();
  EXPECT_EQ(test_output_201_012,ref_output_201_012) << "Failed sparsification with 012->201 nesting order change and 012->012 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(dense_segs, lam021);
  std::string ref_output_201_021 = stream.str();
  auto sparsified_201_021 = sparsify_iteration_space<POL201>(dense_segs, A, perm021);
  stream.clear();
  stream.str("");
  kernel<POL201>(sparsified_201_021, lam_prt);
  std::string test_output_201_021 = stream.str();
  EXPECT_EQ(test_output_201_021,ref_output_201_021) << "Failed sparsification with 012->201 nesting order change and 012->021 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(dense_segs, lam102);
  std::string ref_output_201_102 = stream.str();
  auto sparsified_201_102 = sparsify_iteration_space<POL201>(dense_segs, A, perm102);
  stream.clear();
  stream.str("");
  kernel<POL201>(sparsified_201_102, lam_prt);
  std::string test_output_201_102 = stream.str();
  EXPECT_EQ(test_output_201_102,ref_output_201_102) << "Failed sparsification with 012->201 nesting order change and 012->102 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(dense_segs, lam120);
  std::string ref_output_201_120 = stream.str();
  auto sparsified_201_120 = sparsify_iteration_space<POL201>(dense_segs, A, perm120);
  stream.clear();
  stream.str("");
  kernel<POL201>(sparsified_201_120, lam_prt);
  std::string test_output_201_120 = stream.str();
  EXPECT_EQ(test_output_201_120,ref_output_201_120) << "Failed sparsification with 012->201 nesting order change and 012->120 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(dense_segs, lam201);
  std::string ref_output_201_201 = stream.str();
  auto sparsified_201_201 = sparsify_iteration_space<POL201>(dense_segs, A, perm201);
  stream.clear();
  stream.str("");
  kernel<POL201>(sparsified_201_201, lam_prt);
  std::string test_output_201_201 = stream.str();
  EXPECT_EQ(test_output_201_201,ref_output_201_201) << "Failed sparsification with 012->201 nesting order change and 012->201 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(dense_segs, lam210);
  std::string ref_output_201_210 = stream.str();
  auto sparsified_201_210 = sparsify_iteration_space<POL201>(dense_segs, A, perm210);
  stream.clear();
  stream.str("");
  kernel<POL201>(sparsified_201_210, lam_prt);
  std::string test_output_201_210 = stream.str();
  EXPECT_EQ(test_output_201_210,ref_output_201_210) << "Failed sparsification with 012->201 nesting order change and 012->210 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL210>(dense_segs, lam012);
  std::string ref_output_210_012 = stream.str();
  auto sparsified_210_012 = sparsify_iteration_space<POL210>(dense_segs, A, perm012);
  stream.clear();
  stream.str("");
  kernel<POL210>(sparsified_210_012, lam_prt);
  std::string test_output_210_012 = stream.str();
  EXPECT_EQ(test_output_210_012,ref_output_210_012) << "Failed sparsification with 012->210 nesting order change and 012->012 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL210>(dense_segs, lam021);
  std::string ref_output_210_021 = stream.str();
  auto sparsified_210_021 = sparsify_iteration_space<POL210>(dense_segs, A, perm021);
  stream.clear();
  stream.str("");
  kernel<POL210>(sparsified_210_021, lam_prt);
  std::string test_output_210_021 = stream.str();
  EXPECT_EQ(test_output_210_021,ref_output_210_021) << "Failed sparsification with 012->210 nesting order change and 012->021 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL210>(dense_segs, lam102);
  std::string ref_output_210_102 = stream.str();
  auto sparsified_210_102 = sparsify_iteration_space<POL210>(dense_segs, A, perm102);
  stream.clear();
  stream.str("");
  kernel<POL210>(sparsified_210_102, lam_prt);
  std::string test_output_210_102 = stream.str();
  EXPECT_EQ(test_output_210_102,ref_output_210_102) << "Failed sparsification with 012->210 nesting order change and 012->102 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL210>(dense_segs, lam120);
  std::string ref_output_210_120 = stream.str();
  auto sparsified_210_120 = sparsify_iteration_space<POL210>(dense_segs, A, perm120);
  stream.clear();
  stream.str("");
  kernel<POL210>(sparsified_210_120, lam_prt);
  std::string test_output_210_120 = stream.str();
  EXPECT_EQ(test_output_210_120,ref_output_210_120) << "Failed sparsification with 012->210 nesting order change and 012->120 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL210>(dense_segs, lam201);
  std::string ref_output_210_201 = stream.str();
  auto sparsified_210_201 = sparsify_iteration_space<POL210>(dense_segs, A, perm201);
  stream.clear();
  stream.str("");
  kernel<POL210>(sparsified_210_201, lam_prt);
  std::string test_output_210_201 = stream.str();
  EXPECT_EQ(test_output_210_201,ref_output_210_201) << "Failed sparsification with 012->210 nesting order change and 012->201 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

  A.impl->set_permutation(perm012);
  stream.clear();
  stream.str(std::string());
  kernel<POL210>(dense_segs, lam210);
  std::string ref_output_210_210 = stream.str();
  auto sparsified_210_210 = sparsify_iteration_space<POL210>(dense_segs, A, perm210);
  stream.clear();
  stream.str("");
  kernel<POL210>(sparsified_210_210, lam_prt);
  std::string test_output_210_210 = stream.str();
  EXPECT_EQ(test_output_210_210,ref_output_210_210) << "Failed sparsification with 012->210 nesting order change and 012->210 argument order change\n" << "012 012 reference output: " << ref_output_012_012 << "\n";

 

  auto lam2_01 = [&](auto i0, auto i1, auto i2) {
    if(A2(i0,i1) != 0.0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam2_02 = [&](auto i0, auto i1, auto i2) {
    if(A2(i0,i2) != 0.0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam2_10 = [&](auto i0, auto i1, auto i2) {
    if(A2(i1,i0) != 0.0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam2_21 = [&](auto i0, auto i1, auto i2) {
    if(A2(i2,i1) != 0.0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam2_12 = [&](auto i0, auto i1, auto i2) {
    if(A2(i1,i2) != 0.0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };


  std::vector<camp::idx_t> perm01 = {{0,1}};
  std::vector<camp::idx_t> perm10 = {{1,0}};
  std::vector<camp::idx_t> perm02 = {{0,2}};
  std::vector<camp::idx_t> perm21 = {{2,1}};
  std::vector<camp::idx_t> perm12 = {{1,2}};
/*

  A2.impl->set_permutation(perm01);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam2_01);
  std::string ref_output_2_012_01 = stream.str();
  auto sparsified_2_012_01 = sparsify_iteration_space<POL012>(dense_segs, A2, perm01);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(sparsified_2_012_01, lam_prt);
  std::string test_output_2_012_01 = stream.str();
  EXPECT_EQ(test_output_2_012_01, ref_output_2_012_01) << "Failed sparsification with 012 nesting order and 01 access order for 2D data.\n";

  A2.impl->set_permutation(perm01);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam2_10);
  std::string ref_output_2_012_10 = stream.str();
  auto sparsified_2_012_10 = sparsify_iteration_space<POL012>(dense_segs, A2, perm10);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(sparsified_2_012_10, lam_prt);
  std::string test_output_2_012_10 = stream.str();
  EXPECT_EQ(test_output_2_012_10, ref_output_2_012_10) << "Failed sparsification with 012 nesting order and 10 access order for 2D data.\n";

  A2.impl->set_permutation(perm01);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam2_21);
  std::string ref_output_2_012_21 = stream.str();
  auto sparsified_2_012_21 = sparsify_iteration_space<POL012>(dense_segs, A2, perm21);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(sparsified_2_012_21, lam_prt);
  std::string test_output_2_012_21 = stream.str();
  EXPECT_EQ(test_output_2_012_21, ref_output_2_012_21) << "Failed sparsification with 012 nesting order and 21 access order for 2D data.\n";

  A2.impl->set_permutation(perm01);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(dense_segs, lam2_01);
  std::string ref_output_2_201_01 = stream.str();
  auto sparsified_2_201_01 = sparsify_iteration_space<POL201>(dense_segs, A2, perm01);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(sparsified_2_201_01, lam_prt);
  std::string test_output_2_201_01 = stream.str();
  EXPECT_EQ(test_output_2_201_01, ref_output_2_201_01) << "Failed sparsification with 201 nesting order and 01 access order for 2D data.\n";

  A2.impl->set_permutation(perm01);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(dense_segs, lam2_02);
  std::string ref_output_2_201_02 = stream.str();
  auto sparsified_2_201_02 = sparsify_iteration_space<POL201>(dense_segs, A2, perm02);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(sparsified_2_201_02, lam_prt);
  std::string test_output_2_201_02 = stream.str();
  EXPECT_EQ(test_output_2_201_02, ref_output_2_201_02) << "Failed sparsification with 201 nesting order and 02 access order for 2D data.\n";

  A2.impl->set_permutation(perm01);
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(dense_segs, lam2_21);
  std::string ref_output_2_102_21 = stream.str();
  auto sparsified_2_102_21 = sparsify_iteration_space<POL102>(dense_segs, A2, perm21);
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(sparsified_2_102_21, lam_prt);
  std::string test_output_2_102_21 = stream.str();
  EXPECT_EQ(test_output_2_102_21, ref_output_2_102_21) << "Failed sparsification with 102 nesting order and 21 access order for 2D data.\n";

  A2.impl->set_permutation(perm01);
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(dense_segs, lam2_12);
  std::string ref_output_2_102_12 = stream.str();
  auto sparsified_2_102_12 = sparsify_iteration_space<POL102>(dense_segs, A2, perm12);
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(sparsified_2_102_12, lam_prt);
  std::string test_output_2_102_12 = stream.str();
  EXPECT_EQ(test_output_2_102_12, ref_output_2_102_12) << "Failed sparsification with 102 nesting order and 12 access order for 2D data.\n";
*/
/*
  auto lam4_0012 = [&](auto i0, auto i1, auto i2) {
    if(A4(i0,i0,i1,i2) != 0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };
  auto lam4_2021 = [&](auto i0, auto i1, auto i2) {
    if(A4(i2,i0,i2,i1) != 0) {stream << i0 << "," << i1 << "," << i2 << "  ";}
  };

  std::vector<camp::idx_t> perm0012 = {{0,0,1,2}};
  std::vector<camp::idx_t> perm0123 = {{0,1,2,3}};
  std::vector<camp::idx_t> perm2021 = {{2,0,2,1}};

  A4.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam4_0012);
  std::string ref_output_4_012_0012 = stream.str();
  auto sparsified_4_012_0012 = sparsify_iteration_space<POL012>(dense_segs, A4, perm0012);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(sparsified_4_012_0012, lam_prt);
  std::string test_output_4_012_0012 = stream.str();
  EXPECT_EQ(test_output_4_012_0012, ref_output_4_012_0012) << "Failed sparsification with 012 nesting order and 0012 access order for 4D data.\n";

  A4.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam4_2021);
  std::string ref_output_4_012_2021 = stream.str();
  auto sparsified_4_012_2021 = sparsify_iteration_space<POL012>(dense_segs, A4, perm2021);
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(sparsified_4_012_2021, lam_prt);
  std::string test_output_4_012_2021 = stream.str();
  EXPECT_EQ(test_output_4_012_2021, ref_output_4_012_2021) << "Failed sparsification with 012 nesting order and 2021 access order for 4D data.\n";

  A4.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(dense_segs, lam4_0012);
  std::string ref_output_4_201_0012 = stream.str();
  auto sparsified_4_201_0012 = sparsify_iteration_space<POL201>(dense_segs, A4, perm0012);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(sparsified_4_201_0012, lam_prt);
  std::string test_output_4_201_0012 = stream.str();
  EXPECT_EQ(test_output_4_201_0012, ref_output_4_201_0012) << "Failed sparsification with 201 nesting order and 0012 access order for 4D data.\n";

  A4.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(dense_segs, lam4_2021);
  std::string ref_output_4_201_2021 = stream.str();
  auto sparsified_4_201_2021 = sparsify_iteration_space<POL201>(dense_segs, A4, perm2021);
  stream.clear();
  stream.str(std::string());
  kernel<POL201>(sparsified_4_201_2021, lam_prt);
  std::string test_output_4_201_2021 = stream.str();
  EXPECT_EQ(test_output_4_201_2021, ref_output_4_201_2021) << "Failed sparsification with 201 nesting order and 2021 access order for 4D data.\n";
*/
}


TYPED_TEST_SUITE_P(SparseExampleTest);
template <typename T>
class SparseExampleTest : public ::testing::Test
{
};

TYPED_TEST_P(SparseExampleTest, SparseExample)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  SparseExampleTestSparsify3<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SparseExampleTest, SparseExample);
