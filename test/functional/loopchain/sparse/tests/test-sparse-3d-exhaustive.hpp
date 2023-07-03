#include "RAJA_test-base.hpp"
#include <string>
#include <vector>


template <typename EXEC_POLICY>
void SparseExampleTestSparsify3() {
  using namespace RAJA;
  
  camp::idx_t N = 5;
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

  
  A.impl.set_permutation(perm012);
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
  
  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

  A.impl.set_permutation(perm012);
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

}


TYPED_TEST_SUITE_P(Sparsify3DExhaustiveTest);
template <typename T>
class Sparsify3DExhaustiveTest : public ::testing::Test
{
};

TYPED_TEST_P(Sparsify3DExhaustiveTest, Sparsify3DExhaustive)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  SparseExampleTestSparsify3<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(Sparsify3DExhaustiveTest, Sparsify3DExhaustive);
