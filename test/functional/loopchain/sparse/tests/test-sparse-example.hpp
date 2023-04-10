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
  auto sparsified_012_012 = sparsify_iteration_space_equal_dims<POL012>(dense_segs, A, perm012);
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
  auto sparsified_012_021 = sparsify_iteration_space_equal_dims<POL012>(dense_segs, A, perm021);
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
  auto sparsified_012_102 = sparsify_iteration_space_equal_dims<POL012>(dense_segs, A, perm102);
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
  auto sparsified_012_120 = sparsify_iteration_space_equal_dims<POL012>(dense_segs, A, perm120);
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
  auto sparsified_012_201 = sparsify_iteration_space_equal_dims<POL012>(dense_segs, A, perm201);
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
  auto sparsified_012_210 = sparsify_iteration_space_equal_dims<POL012>(dense_segs, A, perm210);
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
  auto sparsified_021_012 = sparsify_iteration_space_equal_dims<POL021>(dense_segs, A, perm012);
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
  auto sparsified_021_021 = sparsify_iteration_space_equal_dims<POL021>(dense_segs, A, perm021);
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
  auto sparsified_021_102 = sparsify_iteration_space_equal_dims<POL021>(dense_segs, A, perm102);
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
  auto sparsified_021_120 = sparsify_iteration_space_equal_dims<POL021>(dense_segs, A, perm120);
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
  auto sparsified_021_201 = sparsify_iteration_space_equal_dims<POL021>(dense_segs, A, perm201);
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
  auto sparsified_021_210 = sparsify_iteration_space_equal_dims<POL021>(dense_segs, A, perm210);
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
  auto sparsified_102_012 = sparsify_iteration_space_equal_dims<POL102>(dense_segs, A, perm012);
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
  auto sparsified_102_021 = sparsify_iteration_space_equal_dims<POL102>(dense_segs, A, perm021);
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
  auto sparsified_102_102 = sparsify_iteration_space_equal_dims<POL102>(dense_segs, A, perm102);
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
  auto sparsified_102_120 = sparsify_iteration_space_equal_dims<POL102>(dense_segs, A, perm120);
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
  auto sparsified_102_201 = sparsify_iteration_space_equal_dims<POL102>(dense_segs, A, perm201);
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
  auto sparsified_102_210 = sparsify_iteration_space_equal_dims<POL102>(dense_segs, A, perm210);
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
  auto sparsified_120_012 = sparsify_iteration_space_equal_dims<POL120>(dense_segs, A, perm012);
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
  auto sparsified_120_021 = sparsify_iteration_space_equal_dims<POL120>(dense_segs, A, perm021);
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
  auto sparsified_120_102 = sparsify_iteration_space_equal_dims<POL120>(dense_segs, A, perm102);
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
  auto sparsified_120_120 = sparsify_iteration_space_equal_dims<POL120>(dense_segs, A, perm120);
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
  auto sparsified_120_201 = sparsify_iteration_space_equal_dims<POL120>(dense_segs, A, perm201);
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
  auto sparsified_120_210 = sparsify_iteration_space_equal_dims<POL120>(dense_segs, A, perm210);
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
  auto sparsified_201_012 = sparsify_iteration_space_equal_dims<POL201>(dense_segs, A, perm012);
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
  auto sparsified_201_021 = sparsify_iteration_space_equal_dims<POL201>(dense_segs, A, perm021);
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
  auto sparsified_201_102 = sparsify_iteration_space_equal_dims<POL201>(dense_segs, A, perm102);
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
  auto sparsified_201_120 = sparsify_iteration_space_equal_dims<POL201>(dense_segs, A, perm120);
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
  auto sparsified_201_201 = sparsify_iteration_space_equal_dims<POL201>(dense_segs, A, perm201);
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
  auto sparsified_201_210 = sparsify_iteration_space_equal_dims<POL201>(dense_segs, A, perm210);
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
  auto sparsified_210_012 = sparsify_iteration_space_equal_dims<POL210>(dense_segs, A, perm012);
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
  auto sparsified_210_021 = sparsify_iteration_space_equal_dims<POL210>(dense_segs, A, perm021);
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
  auto sparsified_210_102 = sparsify_iteration_space_equal_dims<POL210>(dense_segs, A, perm102);
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
  auto sparsified_210_120 = sparsify_iteration_space_equal_dims<POL210>(dense_segs, A, perm120);
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
  auto sparsified_210_201 = sparsify_iteration_space_equal_dims<POL210>(dense_segs, A, perm201);
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
  auto sparsified_210_210 = sparsify_iteration_space_equal_dims<POL210>(dense_segs, A, perm210);
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

  std::vector<camp::idx_t> perm01 = {{0,1}};
  std::vector<camp::idx_t> perm10 = {{1,0}};
  std::vector<camp::idx_t> perm02 = {{0,2}};
  std::vector<camp::idx_t> perm21 = {{2,1}};


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

}


template <typename EXEC_POLICY>
void SparseExampleTestSparsify4() {
  using namespace RAJA;

  camp::idx_t N = 8;
  auto A = make_random_sparse_view4<double>(N);
  
  auto dense_segs = tuple_repeat<4>(RangeSegment(0,N));

  std::stringstream stream;

  auto lam_prt = [&](auto i0,auto i1,auto i2,auto i3) { stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";};
  auto lam0123 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i0,i1,i2,i3) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam0132 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i0,i1,i3,i2) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam0213 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i0,i2,i1,i3) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam0231 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i0,i2,i3,i1) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam0321 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i0,i3,i2,i1) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam1023 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i1,i0,i2,i3) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam1032 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i1,i0,i3,i2) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam1203 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i1,i2,i0,i3) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam1230 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i1,i2,i3,i0) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam1302 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i1,i3,i0,i2) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam1320 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i1,i3,i2,i0) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam2013 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i2,i0,i1,i3) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam2130 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i2,i1,i3,i0) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam2301 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i2,i3,i0,i1) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam2310 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i2,i3,i1,i0) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam3012 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i3,i0,i1,i2) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam3102 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i3,i1,i0,i2) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam3120 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i3,i1,i2,i0) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam3201 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i3,i2,i0,i1) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  auto lam3210 = [&](auto i0,auto i1,auto i2,auto i3) {if (A(i3,i2,i1,i0) != 0) {stream << i0 << "," << i1 << "," << i2 << "," << i3 << "  ";}};
  using POL0123 = KernelPolicy<statement::For<0,loop_exec,statement::For<1,loop_exec,statement::For<2,loop_exec,statement::For<3,loop_exec,statement::Lambda<0>>>>>>;
  using POL0132 = KernelPolicy<statement::For<0,loop_exec,statement::For<1,loop_exec,statement::For<3,loop_exec,statement::For<2,loop_exec,statement::Lambda<0>>>>>>;
  using POL0231 = KernelPolicy<statement::For<0,loop_exec,statement::For<2,loop_exec,statement::For<3,loop_exec,statement::For<1,loop_exec,statement::Lambda<0>>>>>>;
  using POL0312 = KernelPolicy<statement::For<0,loop_exec,statement::For<3,loop_exec,statement::For<1,loop_exec,statement::For<2,loop_exec,statement::Lambda<0>>>>>>;
  using POL0321 = KernelPolicy<statement::For<0,loop_exec,statement::For<3,loop_exec,statement::For<2,loop_exec,statement::For<1,loop_exec,statement::Lambda<0>>>>>>;
  using POL1023 = KernelPolicy<statement::For<1,loop_exec,statement::For<0,loop_exec,statement::For<2,loop_exec,statement::For<3,loop_exec,statement::Lambda<0>>>>>>;
  using POL1032 = KernelPolicy<statement::For<1,loop_exec,statement::For<0,loop_exec,statement::For<3,loop_exec,statement::For<2,loop_exec,statement::Lambda<0>>>>>>;
  using POL1203 = KernelPolicy<statement::For<1,loop_exec,statement::For<2,loop_exec,statement::For<0,loop_exec,statement::For<3,loop_exec,statement::Lambda<0>>>>>>;
  using POL1230 = KernelPolicy<statement::For<1,loop_exec,statement::For<2,loop_exec,statement::For<3,loop_exec,statement::For<0,loop_exec,statement::Lambda<0>>>>>>;
  using POL1320 = KernelPolicy<statement::For<1,loop_exec,statement::For<3,loop_exec,statement::For<2,loop_exec,statement::For<0,loop_exec,statement::Lambda<0>>>>>>;
  using POL2013 = KernelPolicy<statement::For<2,loop_exec,statement::For<0,loop_exec,statement::For<1,loop_exec,statement::For<3,loop_exec,statement::Lambda<0>>>>>>;
  using POL2031 = KernelPolicy<statement::For<2,loop_exec,statement::For<0,loop_exec,statement::For<3,loop_exec,statement::For<1,loop_exec,statement::Lambda<0>>>>>>;
  using POL2103 = KernelPolicy<statement::For<2,loop_exec,statement::For<1,loop_exec,statement::For<0,loop_exec,statement::For<3,loop_exec,statement::Lambda<0>>>>>>;
  using POL2130 = KernelPolicy<statement::For<2,loop_exec,statement::For<1,loop_exec,statement::For<3,loop_exec,statement::For<0,loop_exec,statement::Lambda<0>>>>>>;
  using POL2310 = KernelPolicy<statement::For<2,loop_exec,statement::For<3,loop_exec,statement::For<1,loop_exec,statement::For<0,loop_exec,statement::Lambda<0>>>>>>;
  using POL3012 = KernelPolicy<statement::For<3,loop_exec,statement::For<0,loop_exec,statement::For<1,loop_exec,statement::For<2,loop_exec,statement::Lambda<0>>>>>>;
  using POL3021 = KernelPolicy<statement::For<3,loop_exec,statement::For<0,loop_exec,statement::For<2,loop_exec,statement::For<1,loop_exec,statement::Lambda<0>>>>>>;
  using POL3102 = KernelPolicy<statement::For<3,loop_exec,statement::For<1,loop_exec,statement::For<0,loop_exec,statement::For<2,loop_exec,statement::Lambda<0>>>>>>;
  using POL3120 = KernelPolicy<statement::For<3,loop_exec,statement::For<1,loop_exec,statement::For<2,loop_exec,statement::For<0,loop_exec,statement::Lambda<0>>>>>>;
  using POL3201 = KernelPolicy<statement::For<3,loop_exec,statement::For<2,loop_exec,statement::For<0,loop_exec,statement::For<1,loop_exec,statement::Lambda<0>>>>>>;
  using POL3210 = KernelPolicy<statement::For<3,loop_exec,statement::For<2,loop_exec,statement::For<1,loop_exec,statement::For<0,loop_exec,statement::Lambda<0>>>>>>;
  std::vector<camp::idx_t> perm0123 = {0,1,2,3};
  std::vector<camp::idx_t> perm0132 = {0,1,3,2};
  std::vector<camp::idx_t> perm0213 = {0,2,1,3};
  std::vector<camp::idx_t> perm0231 = {0,2,3,1};
  std::vector<camp::idx_t> perm0312 = {0,3,1,2};
  std::vector<camp::idx_t> perm0321 = {0,3,2,1};
  std::vector<camp::idx_t> perm1023 = {1,0,2,3};
  std::vector<camp::idx_t> perm1032 = {1,0,3,2};
  std::vector<camp::idx_t> perm1203 = {1,2,0,3};
  std::vector<camp::idx_t> perm1230 = {1,2,3,0};
  std::vector<camp::idx_t> perm1302 = {1,3,0,2};
  std::vector<camp::idx_t> perm1320 = {1,3,2,0};
  std::vector<camp::idx_t> perm2013 = {2,0,1,3};
  std::vector<camp::idx_t> perm2031 = {2,0,3,1};
  std::vector<camp::idx_t> perm2103 = {2,1,0,3};
  std::vector<camp::idx_t> perm2130 = {2,1,3,0};
  std::vector<camp::idx_t> perm2301 = {2,3,0,1};
  std::vector<camp::idx_t> perm2310 = {2,3,1,0};
  std::vector<camp::idx_t> perm3012 = {3,0,1,2};
  std::vector<camp::idx_t> perm3021 = {3,0,2,1};
  std::vector<camp::idx_t> perm3102 = {3,1,0,2};
  std::vector<camp::idx_t> perm3120 = {3,1,2,0};
  std::vector<camp::idx_t> perm3201 = {3,2,0,1};
  std::vector<camp::idx_t> perm3210 = {3,2,1,0};
  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0123>(dense_segs, lam0123);
  std::string ref_output_0123_0123 = stream.str();
  auto sparsified_0123_0123 = sparsify_iteration_space_equal_dims<POL0123>(dense_segs, A, perm0123);
  stream.clear();
  stream.str("");
  kernel<POL0123>(sparsified_0123_0123, lam_prt);
  std::string test_output_0123_0123 = stream.str();
  EXPECT_EQ(test_output_0123_0123,ref_output_0123_0123) << "Failed sparsification with 0123->0123 nesting order change and 0123->0123 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0132>(dense_segs, lam3201);
  std::string ref_output_0132_3201 = stream.str();
  auto sparsified_0132_3201 = sparsify_iteration_space_equal_dims<POL0132>(dense_segs, A, perm3201);
  stream.clear();
  stream.str("");
  kernel<POL0132>(sparsified_0132_3201, lam_prt);
  std::string test_output_0132_3201 = stream.str();
  EXPECT_EQ(test_output_0132_3201,ref_output_0132_3201) << "Failed sparsification with 0123->0132 nesting order change and 0123->3201 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0231>(dense_segs, lam0213);
  std::string ref_output_0231_0213 = stream.str();
  auto sparsified_0231_0213 = sparsify_iteration_space_equal_dims<POL0231>(dense_segs, A, perm0213);
  stream.clear();
  stream.str("");
  kernel<POL0231>(sparsified_0231_0213, lam_prt);
  std::string test_output_0231_0213 = stream.str();
  EXPECT_EQ(test_output_0231_0213,ref_output_0231_0213) << "Failed sparsification with 0123->0231 nesting order change and 0123->0213 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0312>(dense_segs, lam0123);
  std::string ref_output_0312_0123 = stream.str();
  auto sparsified_0312_0123 = sparsify_iteration_space_equal_dims<POL0312>(dense_segs, A, perm0123);
  stream.clear();
  stream.str("");
  kernel<POL0312>(sparsified_0312_0123, lam_prt);
  std::string test_output_0312_0123 = stream.str();
  EXPECT_EQ(test_output_0312_0123,ref_output_0312_0123) << "Failed sparsification with 0123->0312 nesting order change and 0123->0123 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0312>(dense_segs, lam0132);
  std::string ref_output_0312_0132 = stream.str();
  auto sparsified_0312_0132 = sparsify_iteration_space_equal_dims<POL0312>(dense_segs, A, perm0132);
  stream.clear();
  stream.str("");
  kernel<POL0312>(sparsified_0312_0132, lam_prt);
  std::string test_output_0312_0132 = stream.str();
  EXPECT_EQ(test_output_0312_0132,ref_output_0312_0132) << "Failed sparsification with 0123->0312 nesting order change and 0123->0132 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0312>(dense_segs, lam1203);
  std::string ref_output_0312_1203 = stream.str();
  auto sparsified_0312_1203 = sparsify_iteration_space_equal_dims<POL0312>(dense_segs, A, perm1203);
  stream.clear();
  stream.str("");
  kernel<POL0312>(sparsified_0312_1203, lam_prt);
  std::string test_output_0312_1203 = stream.str();
  EXPECT_EQ(test_output_0312_1203,ref_output_0312_1203) << "Failed sparsification with 0123->0312 nesting order change and 0123->1203 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0312>(dense_segs, lam2310);
  std::string ref_output_0312_2310 = stream.str();
  auto sparsified_0312_2310 = sparsify_iteration_space_equal_dims<POL0312>(dense_segs, A, perm2310);
  stream.clear();
  stream.str("");
  kernel<POL0312>(sparsified_0312_2310, lam_prt);
  std::string test_output_0312_2310 = stream.str();
  EXPECT_EQ(test_output_0312_2310,ref_output_0312_2310) << "Failed sparsification with 0123->0312 nesting order change and 0123->2310 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0321>(dense_segs, lam0213);
  std::string ref_output_0321_0213 = stream.str();
  auto sparsified_0321_0213 = sparsify_iteration_space_equal_dims<POL0321>(dense_segs, A, perm0213);
  stream.clear();
  stream.str("");
  kernel<POL0321>(sparsified_0321_0213, lam_prt);
  std::string test_output_0321_0213 = stream.str();
  EXPECT_EQ(test_output_0321_0213,ref_output_0321_0213) << "Failed sparsification with 0123->0321 nesting order change and 0123->0213 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1023>(dense_segs, lam0213);
  std::string ref_output_1023_0213 = stream.str();
  auto sparsified_1023_0213 = sparsify_iteration_space_equal_dims<POL1023>(dense_segs, A, perm0213);
  stream.clear();
  stream.str("");
  kernel<POL1023>(sparsified_1023_0213, lam_prt);
  std::string test_output_1023_0213 = stream.str();
  EXPECT_EQ(test_output_1023_0213,ref_output_1023_0213) << "Failed sparsification with 0123->1023 nesting order change and 0123->0213 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1023>(dense_segs, lam3120);
  std::string ref_output_1023_3120 = stream.str();
  auto sparsified_1023_3120 = sparsify_iteration_space_equal_dims<POL1023>(dense_segs, A, perm3120);
  stream.clear();
  stream.str("");
  kernel<POL1023>(sparsified_1023_3120, lam_prt);
  std::string test_output_1023_3120 = stream.str();
  EXPECT_EQ(test_output_1023_3120,ref_output_1023_3120) << "Failed sparsification with 0123->1023 nesting order change and 0123->3120 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1032>(dense_segs, lam0213);
  std::string ref_output_1032_0213 = stream.str();
  auto sparsified_1032_0213 = sparsify_iteration_space_equal_dims<POL1032>(dense_segs, A, perm0213);
  stream.clear();
  stream.str("");
  kernel<POL1032>(sparsified_1032_0213, lam_prt);
  std::string test_output_1032_0213 = stream.str();
  EXPECT_EQ(test_output_1032_0213,ref_output_1032_0213) << "Failed sparsification with 0123->1032 nesting order change and 0123->0213 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1032>(dense_segs, lam0231);
  std::string ref_output_1032_0231 = stream.str();
  auto sparsified_1032_0231 = sparsify_iteration_space_equal_dims<POL1032>(dense_segs, A, perm0231);
  stream.clear();
  stream.str("");
  kernel<POL1032>(sparsified_1032_0231, lam_prt);
  std::string test_output_1032_0231 = stream.str();
  EXPECT_EQ(test_output_1032_0231,ref_output_1032_0231) << "Failed sparsification with 0123->1032 nesting order change and 0123->0231 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1032>(dense_segs, lam1023);
  std::string ref_output_1032_1023 = stream.str();
  auto sparsified_1032_1023 = sparsify_iteration_space_equal_dims<POL1032>(dense_segs, A, perm1023);
  stream.clear();
  stream.str("");
  kernel<POL1032>(sparsified_1032_1023, lam_prt);
  std::string test_output_1032_1023 = stream.str();
  EXPECT_EQ(test_output_1032_1023,ref_output_1032_1023) << "Failed sparsification with 0123->1032 nesting order change and 0123->1023 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1032>(dense_segs, lam3102);
  std::string ref_output_1032_3102 = stream.str();
  auto sparsified_1032_3102 = sparsify_iteration_space_equal_dims<POL1032>(dense_segs, A, perm3102);
  stream.clear();
  stream.str("");
  kernel<POL1032>(sparsified_1032_3102, lam_prt);
  std::string test_output_1032_3102 = stream.str();
  EXPECT_EQ(test_output_1032_3102,ref_output_1032_3102) << "Failed sparsification with 0123->1032 nesting order change and 0123->3102 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1203>(dense_segs, lam3102);
  std::string ref_output_1203_3102 = stream.str();
  auto sparsified_1203_3102 = sparsify_iteration_space_equal_dims<POL1203>(dense_segs, A, perm3102);
  stream.clear();
  stream.str("");
  kernel<POL1203>(sparsified_1203_3102, lam_prt);
  std::string test_output_1203_3102 = stream.str();
  EXPECT_EQ(test_output_1203_3102,ref_output_1203_3102) << "Failed sparsification with 0123->1203 nesting order change and 0123->3102 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1230>(dense_segs, lam1230);
  std::string ref_output_1230_1230 = stream.str();
  auto sparsified_1230_1230 = sparsify_iteration_space_equal_dims<POL1230>(dense_segs, A, perm1230);
  stream.clear();
  stream.str("");
  kernel<POL1230>(sparsified_1230_1230, lam_prt);
  std::string test_output_1230_1230 = stream.str();
  EXPECT_EQ(test_output_1230_1230,ref_output_1230_1230) << "Failed sparsification with 0123->1230 nesting order change and 0123->1230 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1320>(dense_segs, lam2013);
  std::string ref_output_1320_2013 = stream.str();
  auto sparsified_1320_2013 = sparsify_iteration_space_equal_dims<POL1320>(dense_segs, A, perm2013);
  stream.clear();
  stream.str("");
  kernel<POL1320>(sparsified_1320_2013, lam_prt);
  std::string test_output_1320_2013 = stream.str();
  EXPECT_EQ(test_output_1320_2013,ref_output_1320_2013) << "Failed sparsification with 0123->1320 nesting order change and 0123->2013 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2013>(dense_segs, lam0123);
  std::string ref_output_2013_0123 = stream.str();
  auto sparsified_2013_0123 = sparsify_iteration_space_equal_dims<POL2013>(dense_segs, A, perm0123);
  stream.clear();
  stream.str("");
  kernel<POL2013>(sparsified_2013_0123, lam_prt);
  std::string test_output_2013_0123 = stream.str();
  EXPECT_EQ(test_output_2013_0123,ref_output_2013_0123) << "Failed sparsification with 0123->2013 nesting order change and 0123->0123 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2013>(dense_segs, lam1320);
  std::string ref_output_2013_1320 = stream.str();
  auto sparsified_2013_1320 = sparsify_iteration_space_equal_dims<POL2013>(dense_segs, A, perm1320);
  stream.clear();
  stream.str("");
  kernel<POL2013>(sparsified_2013_1320, lam_prt);
  std::string test_output_2013_1320 = stream.str();
  EXPECT_EQ(test_output_2013_1320,ref_output_2013_1320) << "Failed sparsification with 0123->2013 nesting order change and 0123->1320 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2031>(dense_segs, lam0321);
  std::string ref_output_2031_0321 = stream.str();
  auto sparsified_2031_0321 = sparsify_iteration_space_equal_dims<POL2031>(dense_segs, A, perm0321);
  stream.clear();
  stream.str("");
  kernel<POL2031>(sparsified_2031_0321, lam_prt);
  std::string test_output_2031_0321 = stream.str();
  EXPECT_EQ(test_output_2031_0321,ref_output_2031_0321) << "Failed sparsification with 0123->2031 nesting order change and 0123->0321 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2031>(dense_segs, lam3012);
  std::string ref_output_2031_3012 = stream.str();
  auto sparsified_2031_3012 = sparsify_iteration_space_equal_dims<POL2031>(dense_segs, A, perm3012);
  stream.clear();
  stream.str("");
  kernel<POL2031>(sparsified_2031_3012, lam_prt);
  std::string test_output_2031_3012 = stream.str();
  EXPECT_EQ(test_output_2031_3012,ref_output_2031_3012) << "Failed sparsification with 0123->2031 nesting order change and 0123->3012 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2103>(dense_segs, lam0123);
  std::string ref_output_2103_0123 = stream.str();
  auto sparsified_2103_0123 = sparsify_iteration_space_equal_dims<POL2103>(dense_segs, A, perm0123);
  stream.clear();
  stream.str("");
  kernel<POL2103>(sparsified_2103_0123, lam_prt);
  std::string test_output_2103_0123 = stream.str();
  EXPECT_EQ(test_output_2103_0123,ref_output_2103_0123) << "Failed sparsification with 0123->2103 nesting order change and 0123->0123 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2103>(dense_segs, lam3012);
  std::string ref_output_2103_3012 = stream.str();
  auto sparsified_2103_3012 = sparsify_iteration_space_equal_dims<POL2103>(dense_segs, A, perm3012);
  stream.clear();
  stream.str("");
  kernel<POL2103>(sparsified_2103_3012, lam_prt);
  std::string test_output_2103_3012 = stream.str();
  EXPECT_EQ(test_output_2103_3012,ref_output_2103_3012) << "Failed sparsification with 0123->2103 nesting order change and 0123->3012 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2130>(dense_segs, lam2013);
  std::string ref_output_2130_2013 = stream.str();
  auto sparsified_2130_2013 = sparsify_iteration_space_equal_dims<POL2130>(dense_segs, A, perm2013);
  stream.clear();
  stream.str("");
  kernel<POL2130>(sparsified_2130_2013, lam_prt);
  std::string test_output_2130_2013 = stream.str();
  EXPECT_EQ(test_output_2130_2013,ref_output_2130_2013) << "Failed sparsification with 0123->2130 nesting order change and 0123->2013 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2130>(dense_segs, lam2310);
  std::string ref_output_2130_2310 = stream.str();
  auto sparsified_2130_2310 = sparsify_iteration_space_equal_dims<POL2130>(dense_segs, A, perm2310);
  stream.clear();
  stream.str("");
  kernel<POL2130>(sparsified_2130_2310, lam_prt);
  std::string test_output_2130_2310 = stream.str();
  EXPECT_EQ(test_output_2130_2310,ref_output_2130_2310) << "Failed sparsification with 0123->2130 nesting order change and 0123->2310 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2310>(dense_segs, lam1203);
  std::string ref_output_2310_1203 = stream.str();
  auto sparsified_2310_1203 = sparsify_iteration_space_equal_dims<POL2310>(dense_segs, A, perm1203);
  stream.clear();
  stream.str("");
  kernel<POL2310>(sparsified_2310_1203, lam_prt);
  std::string test_output_2310_1203 = stream.str();
  EXPECT_EQ(test_output_2310_1203,ref_output_2310_1203) << "Failed sparsification with 0123->2310 nesting order change and 0123->1203 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2310>(dense_segs, lam1320);
  std::string ref_output_2310_1320 = stream.str();
  auto sparsified_2310_1320 = sparsify_iteration_space_equal_dims<POL2310>(dense_segs, A, perm1320);
  stream.clear();
  stream.str("");
  kernel<POL2310>(sparsified_2310_1320, lam_prt);
  std::string test_output_2310_1320 = stream.str();
  EXPECT_EQ(test_output_2310_1320,ref_output_2310_1320) << "Failed sparsification with 0123->2310 nesting order change and 0123->1320 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2310>(dense_segs, lam2130);
  std::string ref_output_2310_2130 = stream.str();
  auto sparsified_2310_2130 = sparsify_iteration_space_equal_dims<POL2310>(dense_segs, A, perm2130);
  stream.clear();
  stream.str("");
  kernel<POL2310>(sparsified_2310_2130, lam_prt);
  std::string test_output_2310_2130 = stream.str();
  EXPECT_EQ(test_output_2310_2130,ref_output_2310_2130) << "Failed sparsification with 0123->2310 nesting order change and 0123->2130 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3012>(dense_segs, lam1032);
  std::string ref_output_3012_1032 = stream.str();
  auto sparsified_3012_1032 = sparsify_iteration_space_equal_dims<POL3012>(dense_segs, A, perm1032);
  stream.clear();
  stream.str("");
  kernel<POL3012>(sparsified_3012_1032, lam_prt);
  std::string test_output_3012_1032 = stream.str();
  EXPECT_EQ(test_output_3012_1032,ref_output_3012_1032) << "Failed sparsification with 0123->3012 nesting order change and 0123->1032 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3012>(dense_segs, lam2130);
  std::string ref_output_3012_2130 = stream.str();
  auto sparsified_3012_2130 = sparsify_iteration_space_equal_dims<POL3012>(dense_segs, A, perm2130);
  stream.clear();
  stream.str("");
  kernel<POL3012>(sparsified_3012_2130, lam_prt);
  std::string test_output_3012_2130 = stream.str();
  EXPECT_EQ(test_output_3012_2130,ref_output_3012_2130) << "Failed sparsification with 0123->3012 nesting order change and 0123->2130 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3021>(dense_segs, lam1203);
  std::string ref_output_3021_1203 = stream.str();
  auto sparsified_3021_1203 = sparsify_iteration_space_equal_dims<POL3021>(dense_segs, A, perm1203);
  stream.clear();
  stream.str("");
  kernel<POL3021>(sparsified_3021_1203, lam_prt);
  std::string test_output_3021_1203 = stream.str();
  EXPECT_EQ(test_output_3021_1203,ref_output_3021_1203) << "Failed sparsification with 0123->3021 nesting order change and 0123->1203 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3021>(dense_segs, lam3120);
  std::string ref_output_3021_3120 = stream.str();
  auto sparsified_3021_3120 = sparsify_iteration_space_equal_dims<POL3021>(dense_segs, A, perm3120);
  stream.clear();
  stream.str("");
  kernel<POL3021>(sparsified_3021_3120, lam_prt);
  std::string test_output_3021_3120 = stream.str();
  EXPECT_EQ(test_output_3021_3120,ref_output_3021_3120) << "Failed sparsification with 0123->3021 nesting order change and 0123->3120 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3102>(dense_segs, lam3102);
  std::string ref_output_3102_3102 = stream.str();
  auto sparsified_3102_3102 = sparsify_iteration_space_equal_dims<POL3102>(dense_segs, A, perm3102);
  stream.clear();
  stream.str("");
  kernel<POL3102>(sparsified_3102_3102, lam_prt);
  std::string test_output_3102_3102 = stream.str();
  EXPECT_EQ(test_output_3102_3102,ref_output_3102_3102) << "Failed sparsification with 0123->3102 nesting order change and 0123->3102 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3120>(dense_segs, lam0213);
  std::string ref_output_3120_0213 = stream.str();
  auto sparsified_3120_0213 = sparsify_iteration_space_equal_dims<POL3120>(dense_segs, A, perm0213);
  stream.clear();
  stream.str("");
  kernel<POL3120>(sparsified_3120_0213, lam_prt);
  std::string test_output_3120_0213 = stream.str();
  EXPECT_EQ(test_output_3120_0213,ref_output_3120_0213) << "Failed sparsification with 0123->3120 nesting order change and 0123->0213 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3120>(dense_segs, lam3120);
  std::string ref_output_3120_3120 = stream.str();
  auto sparsified_3120_3120 = sparsify_iteration_space_equal_dims<POL3120>(dense_segs, A, perm3120);
  stream.clear();
  stream.str("");
  kernel<POL3120>(sparsified_3120_3120, lam_prt);
  std::string test_output_3120_3120 = stream.str();
  EXPECT_EQ(test_output_3120_3120,ref_output_3120_3120) << "Failed sparsification with 0123->3120 nesting order change and 0123->3120 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3120>(dense_segs, lam3210);
  std::string ref_output_3120_3210 = stream.str();
  auto sparsified_3120_3210 = sparsify_iteration_space_equal_dims<POL3120>(dense_segs, A, perm3210);
  stream.clear();
  stream.str("");
  kernel<POL3120>(sparsified_3120_3210, lam_prt);
  std::string test_output_3120_3210 = stream.str();
  EXPECT_EQ(test_output_3120_3210,ref_output_3120_3210) << "Failed sparsification with 0123->3120 nesting order change and 0123->3210 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3201>(dense_segs, lam1302);
  std::string ref_output_3201_1302 = stream.str();
  auto sparsified_3201_1302 = sparsify_iteration_space_equal_dims<POL3201>(dense_segs, A, perm1302);
  stream.clear();
  stream.str("");
  kernel<POL3201>(sparsified_3201_1302, lam_prt);
  std::string test_output_3201_1302 = stream.str();
  EXPECT_EQ(test_output_3201_1302,ref_output_3201_1302) << "Failed sparsification with 0123->3201 nesting order change and 0123->1302 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3201>(dense_segs, lam3201);
  std::string ref_output_3201_3201 = stream.str();
  auto sparsified_3201_3201 = sparsify_iteration_space_equal_dims<POL3201>(dense_segs, A, perm3201);
  stream.clear();
  stream.str("");
  kernel<POL3201>(sparsified_3201_3201, lam_prt);
  std::string test_output_3201_3201 = stream.str();
  EXPECT_EQ(test_output_3201_3201,ref_output_3201_3201) << "Failed sparsification with 0123->3201 nesting order change and 0123->3201 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3210>(dense_segs, lam2301);
  std::string ref_output_3210_2301 = stream.str();
  auto sparsified_3210_2301 = sparsify_iteration_space_equal_dims<POL3210>(dense_segs, A, perm2301);
  stream.clear();
  stream.str("");
  kernel<POL3210>(sparsified_3210_2301, lam_prt);
  std::string test_output_3210_2301 = stream.str();
  EXPECT_EQ(test_output_3210_2301,ref_output_3210_2301) << "Failed sparsification with 0123->3210 nesting order change and 0123->2301 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3210>(dense_segs, lam2310);
  std::string ref_output_3210_2310 = stream.str();
  auto sparsified_3210_2310 = sparsify_iteration_space_equal_dims<POL3210>(dense_segs, A, perm2310);
  stream.clear();
  stream.str("");
  kernel<POL3210>(sparsified_3210_2310, lam_prt);
  std::string test_output_3210_2310 = stream.str();
  EXPECT_EQ(test_output_3210_2310,ref_output_3210_2310) << "Failed sparsification with 0123->3210 nesting order change and 0123->2310 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl->set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3210>(dense_segs, lam3102);
  std::string ref_output_3210_3102 = stream.str();
  auto sparsified_3210_3102 = sparsify_iteration_space_equal_dims<POL3210>(dense_segs, A, perm3102);
  stream.clear();
  stream.str("");
  kernel<POL3210>(sparsified_3210_3102, lam_prt);
  std::string test_output_3210_3102 = stream.str();
  EXPECT_EQ(test_output_3210_3102,ref_output_3210_3102) << "Failed sparsification with 0123->3210 nesting order change and 0123->3102 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";


  
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
  SparseExampleTestSparsify4<EXEC_POLICY>();
  //SparseExampleTestImpl<EXEC_POLICY>();
  //SparseExampleTest3Impl<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SparseExampleTest, SparseExample);
