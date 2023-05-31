#include "RAJA_test-base.hpp"
#include <string>
#include <vector>

template <typename EXEC_POLICY>
void SparseTestSparsify4() {
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
  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0123>(dense_segs, lam0123);
  std::string ref_output_0123_0123 = stream.str();
  auto sparsified_0123_0123 = sparsify_iteration_space<POL0123>(dense_segs, A, perm0123);
  stream.clear();
  stream.str("");
  kernel<POL0123>(sparsified_0123_0123, lam_prt);
  std::string test_output_0123_0123 = stream.str();
  EXPECT_EQ(test_output_0123_0123,ref_output_0123_0123) << "Failed sparsification with 0123->0123 nesting order change and 0123->0123 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0132>(dense_segs, lam3201);
  std::string ref_output_0132_3201 = stream.str();
  auto sparsified_0132_3201 = sparsify_iteration_space<POL0132>(dense_segs, A, perm3201);
  stream.clear();
  stream.str("");
  kernel<POL0132>(sparsified_0132_3201, lam_prt);
  std::string test_output_0132_3201 = stream.str();
  EXPECT_EQ(test_output_0132_3201,ref_output_0132_3201) << "Failed sparsification with 0123->0132 nesting order change and 0123->3201 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0231>(dense_segs, lam0213);
  std::string ref_output_0231_0213 = stream.str();
  auto sparsified_0231_0213 = sparsify_iteration_space<POL0231>(dense_segs, A, perm0213);
  stream.clear();
  stream.str("");
  kernel<POL0231>(sparsified_0231_0213, lam_prt);
  std::string test_output_0231_0213 = stream.str();
  EXPECT_EQ(test_output_0231_0213,ref_output_0231_0213) << "Failed sparsification with 0123->0231 nesting order change and 0123->0213 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0312>(dense_segs, lam0123);
  std::string ref_output_0312_0123 = stream.str();
  auto sparsified_0312_0123 = sparsify_iteration_space<POL0312>(dense_segs, A, perm0123);
  stream.clear();
  stream.str("");
  kernel<POL0312>(sparsified_0312_0123, lam_prt);
  std::string test_output_0312_0123 = stream.str();
  EXPECT_EQ(test_output_0312_0123,ref_output_0312_0123) << "Failed sparsification with 0123->0312 nesting order change and 0123->0123 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0312>(dense_segs, lam0132);
  std::string ref_output_0312_0132 = stream.str();
  auto sparsified_0312_0132 = sparsify_iteration_space<POL0312>(dense_segs, A, perm0132);
  stream.clear();
  stream.str("");
  kernel<POL0312>(sparsified_0312_0132, lam_prt);
  std::string test_output_0312_0132 = stream.str();
  EXPECT_EQ(test_output_0312_0132,ref_output_0312_0132) << "Failed sparsification with 0123->0312 nesting order change and 0123->0132 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0312>(dense_segs, lam1203);
  std::string ref_output_0312_1203 = stream.str();
  auto sparsified_0312_1203 = sparsify_iteration_space<POL0312>(dense_segs, A, perm1203);
  stream.clear();
  stream.str("");
  kernel<POL0312>(sparsified_0312_1203, lam_prt);
  std::string test_output_0312_1203 = stream.str();
  EXPECT_EQ(test_output_0312_1203,ref_output_0312_1203) << "Failed sparsification with 0123->0312 nesting order change and 0123->1203 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0312>(dense_segs, lam2310);
  std::string ref_output_0312_2310 = stream.str();
  auto sparsified_0312_2310 = sparsify_iteration_space<POL0312>(dense_segs, A, perm2310);
  stream.clear();
  stream.str("");
  kernel<POL0312>(sparsified_0312_2310, lam_prt);
  std::string test_output_0312_2310 = stream.str();
  EXPECT_EQ(test_output_0312_2310,ref_output_0312_2310) << "Failed sparsification with 0123->0312 nesting order change and 0123->2310 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL0321>(dense_segs, lam0213);
  std::string ref_output_0321_0213 = stream.str();
  auto sparsified_0321_0213 = sparsify_iteration_space<POL0321>(dense_segs, A, perm0213);
  stream.clear();
  stream.str("");
  kernel<POL0321>(sparsified_0321_0213, lam_prt);
  std::string test_output_0321_0213 = stream.str();
  EXPECT_EQ(test_output_0321_0213,ref_output_0321_0213) << "Failed sparsification with 0123->0321 nesting order change and 0123->0213 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1023>(dense_segs, lam0213);
  std::string ref_output_1023_0213 = stream.str();
  auto sparsified_1023_0213 = sparsify_iteration_space<POL1023>(dense_segs, A, perm0213);
  stream.clear();
  stream.str("");
  kernel<POL1023>(sparsified_1023_0213, lam_prt);
  std::string test_output_1023_0213 = stream.str();
  EXPECT_EQ(test_output_1023_0213,ref_output_1023_0213) << "Failed sparsification with 0123->1023 nesting order change and 0123->0213 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1023>(dense_segs, lam3120);
  std::string ref_output_1023_3120 = stream.str();
  auto sparsified_1023_3120 = sparsify_iteration_space<POL1023>(dense_segs, A, perm3120);
  stream.clear();
  stream.str("");
  kernel<POL1023>(sparsified_1023_3120, lam_prt);
  std::string test_output_1023_3120 = stream.str();
  EXPECT_EQ(test_output_1023_3120,ref_output_1023_3120) << "Failed sparsification with 0123->1023 nesting order change and 0123->3120 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1032>(dense_segs, lam0213);
  std::string ref_output_1032_0213 = stream.str();
  auto sparsified_1032_0213 = sparsify_iteration_space<POL1032>(dense_segs, A, perm0213);
  stream.clear();
  stream.str("");
  kernel<POL1032>(sparsified_1032_0213, lam_prt);
  std::string test_output_1032_0213 = stream.str();
  EXPECT_EQ(test_output_1032_0213,ref_output_1032_0213) << "Failed sparsification with 0123->1032 nesting order change and 0123->0213 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1032>(dense_segs, lam0231);
  std::string ref_output_1032_0231 = stream.str();
  auto sparsified_1032_0231 = sparsify_iteration_space<POL1032>(dense_segs, A, perm0231);
  stream.clear();
  stream.str("");
  kernel<POL1032>(sparsified_1032_0231, lam_prt);
  std::string test_output_1032_0231 = stream.str();
  EXPECT_EQ(test_output_1032_0231,ref_output_1032_0231) << "Failed sparsification with 0123->1032 nesting order change and 0123->0231 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1032>(dense_segs, lam1023);
  std::string ref_output_1032_1023 = stream.str();
  auto sparsified_1032_1023 = sparsify_iteration_space<POL1032>(dense_segs, A, perm1023);
  stream.clear();
  stream.str("");
  kernel<POL1032>(sparsified_1032_1023, lam_prt);
  std::string test_output_1032_1023 = stream.str();
  EXPECT_EQ(test_output_1032_1023,ref_output_1032_1023) << "Failed sparsification with 0123->1032 nesting order change and 0123->1023 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1032>(dense_segs, lam3102);
  std::string ref_output_1032_3102 = stream.str();
  auto sparsified_1032_3102 = sparsify_iteration_space<POL1032>(dense_segs, A, perm3102);
  stream.clear();
  stream.str("");
  kernel<POL1032>(sparsified_1032_3102, lam_prt);
  std::string test_output_1032_3102 = stream.str();
  EXPECT_EQ(test_output_1032_3102,ref_output_1032_3102) << "Failed sparsification with 0123->1032 nesting order change and 0123->3102 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1203>(dense_segs, lam3102);
  std::string ref_output_1203_3102 = stream.str();
  auto sparsified_1203_3102 = sparsify_iteration_space<POL1203>(dense_segs, A, perm3102);
  stream.clear();
  stream.str("");
  kernel<POL1203>(sparsified_1203_3102, lam_prt);
  std::string test_output_1203_3102 = stream.str();
  EXPECT_EQ(test_output_1203_3102,ref_output_1203_3102) << "Failed sparsification with 0123->1203 nesting order change and 0123->3102 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1230>(dense_segs, lam1230);
  std::string ref_output_1230_1230 = stream.str();
  auto sparsified_1230_1230 = sparsify_iteration_space<POL1230>(dense_segs, A, perm1230);
  stream.clear();
  stream.str("");
  kernel<POL1230>(sparsified_1230_1230, lam_prt);
  std::string test_output_1230_1230 = stream.str();
  EXPECT_EQ(test_output_1230_1230,ref_output_1230_1230) << "Failed sparsification with 0123->1230 nesting order change and 0123->1230 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL1320>(dense_segs, lam2013);
  std::string ref_output_1320_2013 = stream.str();
  auto sparsified_1320_2013 = sparsify_iteration_space<POL1320>(dense_segs, A, perm2013);
  stream.clear();
  stream.str("");
  kernel<POL1320>(sparsified_1320_2013, lam_prt);
  std::string test_output_1320_2013 = stream.str();
  EXPECT_EQ(test_output_1320_2013,ref_output_1320_2013) << "Failed sparsification with 0123->1320 nesting order change and 0123->2013 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2013>(dense_segs, lam0123);
  std::string ref_output_2013_0123 = stream.str();
  auto sparsified_2013_0123 = sparsify_iteration_space<POL2013>(dense_segs, A, perm0123);
  stream.clear();
  stream.str("");
  kernel<POL2013>(sparsified_2013_0123, lam_prt);
  std::string test_output_2013_0123 = stream.str();
  EXPECT_EQ(test_output_2013_0123,ref_output_2013_0123) << "Failed sparsification with 0123->2013 nesting order change and 0123->0123 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2013>(dense_segs, lam1320);
  std::string ref_output_2013_1320 = stream.str();
  auto sparsified_2013_1320 = sparsify_iteration_space<POL2013>(dense_segs, A, perm1320);
  stream.clear();
  stream.str("");
  kernel<POL2013>(sparsified_2013_1320, lam_prt);
  std::string test_output_2013_1320 = stream.str();
  EXPECT_EQ(test_output_2013_1320,ref_output_2013_1320) << "Failed sparsification with 0123->2013 nesting order change and 0123->1320 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2031>(dense_segs, lam0321);
  std::string ref_output_2031_0321 = stream.str();
  auto sparsified_2031_0321 = sparsify_iteration_space<POL2031>(dense_segs, A, perm0321);
  stream.clear();
  stream.str("");
  kernel<POL2031>(sparsified_2031_0321, lam_prt);
  std::string test_output_2031_0321 = stream.str();
  EXPECT_EQ(test_output_2031_0321,ref_output_2031_0321) << "Failed sparsification with 0123->2031 nesting order change and 0123->0321 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2031>(dense_segs, lam3012);
  std::string ref_output_2031_3012 = stream.str();
  auto sparsified_2031_3012 = sparsify_iteration_space<POL2031>(dense_segs, A, perm3012);
  stream.clear();
  stream.str("");
  kernel<POL2031>(sparsified_2031_3012, lam_prt);
  std::string test_output_2031_3012 = stream.str();
  EXPECT_EQ(test_output_2031_3012,ref_output_2031_3012) << "Failed sparsification with 0123->2031 nesting order change and 0123->3012 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2103>(dense_segs, lam0123);
  std::string ref_output_2103_0123 = stream.str();
  auto sparsified_2103_0123 = sparsify_iteration_space<POL2103>(dense_segs, A, perm0123);
  stream.clear();
  stream.str("");
  kernel<POL2103>(sparsified_2103_0123, lam_prt);
  std::string test_output_2103_0123 = stream.str();
  EXPECT_EQ(test_output_2103_0123,ref_output_2103_0123) << "Failed sparsification with 0123->2103 nesting order change and 0123->0123 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2103>(dense_segs, lam3012);
  std::string ref_output_2103_3012 = stream.str();
  auto sparsified_2103_3012 = sparsify_iteration_space<POL2103>(dense_segs, A, perm3012);
  stream.clear();
  stream.str("");
  kernel<POL2103>(sparsified_2103_3012, lam_prt);
  std::string test_output_2103_3012 = stream.str();
  EXPECT_EQ(test_output_2103_3012,ref_output_2103_3012) << "Failed sparsification with 0123->2103 nesting order change and 0123->3012 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2130>(dense_segs, lam2013);
  std::string ref_output_2130_2013 = stream.str();
  auto sparsified_2130_2013 = sparsify_iteration_space<POL2130>(dense_segs, A, perm2013);
  stream.clear();
  stream.str("");
  kernel<POL2130>(sparsified_2130_2013, lam_prt);
  std::string test_output_2130_2013 = stream.str();
  EXPECT_EQ(test_output_2130_2013,ref_output_2130_2013) << "Failed sparsification with 0123->2130 nesting order change and 0123->2013 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2130>(dense_segs, lam2310);
  std::string ref_output_2130_2310 = stream.str();
  auto sparsified_2130_2310 = sparsify_iteration_space<POL2130>(dense_segs, A, perm2310);
  stream.clear();
  stream.str("");
  kernel<POL2130>(sparsified_2130_2310, lam_prt);
  std::string test_output_2130_2310 = stream.str();
  EXPECT_EQ(test_output_2130_2310,ref_output_2130_2310) << "Failed sparsification with 0123->2130 nesting order change and 0123->2310 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2310>(dense_segs, lam1203);
  std::string ref_output_2310_1203 = stream.str();
  auto sparsified_2310_1203 = sparsify_iteration_space<POL2310>(dense_segs, A, perm1203);
  stream.clear();
  stream.str("");
  kernel<POL2310>(sparsified_2310_1203, lam_prt);
  std::string test_output_2310_1203 = stream.str();
  EXPECT_EQ(test_output_2310_1203,ref_output_2310_1203) << "Failed sparsification with 0123->2310 nesting order change and 0123->1203 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2310>(dense_segs, lam1320);
  std::string ref_output_2310_1320 = stream.str();
  auto sparsified_2310_1320 = sparsify_iteration_space<POL2310>(dense_segs, A, perm1320);
  stream.clear();
  stream.str("");
  kernel<POL2310>(sparsified_2310_1320, lam_prt);
  std::string test_output_2310_1320 = stream.str();
  EXPECT_EQ(test_output_2310_1320,ref_output_2310_1320) << "Failed sparsification with 0123->2310 nesting order change and 0123->1320 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL2310>(dense_segs, lam2130);
  std::string ref_output_2310_2130 = stream.str();
  auto sparsified_2310_2130 = sparsify_iteration_space<POL2310>(dense_segs, A, perm2130);
  stream.clear();
  stream.str("");
  kernel<POL2310>(sparsified_2310_2130, lam_prt);
  std::string test_output_2310_2130 = stream.str();
  EXPECT_EQ(test_output_2310_2130,ref_output_2310_2130) << "Failed sparsification with 0123->2310 nesting order change and 0123->2130 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3012>(dense_segs, lam1032);
  std::string ref_output_3012_1032 = stream.str();
  auto sparsified_3012_1032 = sparsify_iteration_space<POL3012>(dense_segs, A, perm1032);
  stream.clear();
  stream.str("");
  kernel<POL3012>(sparsified_3012_1032, lam_prt);
  std::string test_output_3012_1032 = stream.str();
  EXPECT_EQ(test_output_3012_1032,ref_output_3012_1032) << "Failed sparsification with 0123->3012 nesting order change and 0123->1032 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3012>(dense_segs, lam2130);
  std::string ref_output_3012_2130 = stream.str();
  auto sparsified_3012_2130 = sparsify_iteration_space<POL3012>(dense_segs, A, perm2130);
  stream.clear();
  stream.str("");
  kernel<POL3012>(sparsified_3012_2130, lam_prt);
  std::string test_output_3012_2130 = stream.str();
  EXPECT_EQ(test_output_3012_2130,ref_output_3012_2130) << "Failed sparsification with 0123->3012 nesting order change and 0123->2130 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3021>(dense_segs, lam1203);
  std::string ref_output_3021_1203 = stream.str();
  auto sparsified_3021_1203 = sparsify_iteration_space<POL3021>(dense_segs, A, perm1203);
  stream.clear();
  stream.str("");
  kernel<POL3021>(sparsified_3021_1203, lam_prt);
  std::string test_output_3021_1203 = stream.str();
  EXPECT_EQ(test_output_3021_1203,ref_output_3021_1203) << "Failed sparsification with 0123->3021 nesting order change and 0123->1203 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3021>(dense_segs, lam3120);
  std::string ref_output_3021_3120 = stream.str();
  auto sparsified_3021_3120 = sparsify_iteration_space<POL3021>(dense_segs, A, perm3120);
  stream.clear();
  stream.str("");
  kernel<POL3021>(sparsified_3021_3120, lam_prt);
  std::string test_output_3021_3120 = stream.str();
  EXPECT_EQ(test_output_3021_3120,ref_output_3021_3120) << "Failed sparsification with 0123->3021 nesting order change and 0123->3120 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3102>(dense_segs, lam3102);
  std::string ref_output_3102_3102 = stream.str();
  auto sparsified_3102_3102 = sparsify_iteration_space<POL3102>(dense_segs, A, perm3102);
  stream.clear();
  stream.str("");
  kernel<POL3102>(sparsified_3102_3102, lam_prt);
  std::string test_output_3102_3102 = stream.str();
  EXPECT_EQ(test_output_3102_3102,ref_output_3102_3102) << "Failed sparsification with 0123->3102 nesting order change and 0123->3102 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3120>(dense_segs, lam0213);
  std::string ref_output_3120_0213 = stream.str();
  auto sparsified_3120_0213 = sparsify_iteration_space<POL3120>(dense_segs, A, perm0213);
  stream.clear();
  stream.str("");
  kernel<POL3120>(sparsified_3120_0213, lam_prt);
  std::string test_output_3120_0213 = stream.str();
  EXPECT_EQ(test_output_3120_0213,ref_output_3120_0213) << "Failed sparsification with 0123->3120 nesting order change and 0123->0213 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3120>(dense_segs, lam3120);
  std::string ref_output_3120_3120 = stream.str();
  auto sparsified_3120_3120 = sparsify_iteration_space<POL3120>(dense_segs, A, perm3120);
  stream.clear();
  stream.str("");
  kernel<POL3120>(sparsified_3120_3120, lam_prt);
  std::string test_output_3120_3120 = stream.str();
  EXPECT_EQ(test_output_3120_3120,ref_output_3120_3120) << "Failed sparsification with 0123->3120 nesting order change and 0123->3120 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3120>(dense_segs, lam3210);
  std::string ref_output_3120_3210 = stream.str();
  auto sparsified_3120_3210 = sparsify_iteration_space<POL3120>(dense_segs, A, perm3210);
  stream.clear();
  stream.str("");
  kernel<POL3120>(sparsified_3120_3210, lam_prt);
  std::string test_output_3120_3210 = stream.str();
  EXPECT_EQ(test_output_3120_3210,ref_output_3120_3210) << "Failed sparsification with 0123->3120 nesting order change and 0123->3210 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3201>(dense_segs, lam1302);
  std::string ref_output_3201_1302 = stream.str();
  auto sparsified_3201_1302 = sparsify_iteration_space<POL3201>(dense_segs, A, perm1302);
  stream.clear();
  stream.str("");
  kernel<POL3201>(sparsified_3201_1302, lam_prt);
  std::string test_output_3201_1302 = stream.str();
  EXPECT_EQ(test_output_3201_1302,ref_output_3201_1302) << "Failed sparsification with 0123->3201 nesting order change and 0123->1302 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3201>(dense_segs, lam3201);
  std::string ref_output_3201_3201 = stream.str();
  auto sparsified_3201_3201 = sparsify_iteration_space<POL3201>(dense_segs, A, perm3201);
  stream.clear();
  stream.str("");
  kernel<POL3201>(sparsified_3201_3201, lam_prt);
  std::string test_output_3201_3201 = stream.str();
  EXPECT_EQ(test_output_3201_3201,ref_output_3201_3201) << "Failed sparsification with 0123->3201 nesting order change and 0123->3201 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3210>(dense_segs, lam2301);
  std::string ref_output_3210_2301 = stream.str();
  auto sparsified_3210_2301 = sparsify_iteration_space<POL3210>(dense_segs, A, perm2301);
  stream.clear();
  stream.str("");
  kernel<POL3210>(sparsified_3210_2301, lam_prt);
  std::string test_output_3210_2301 = stream.str();
  EXPECT_EQ(test_output_3210_2301,ref_output_3210_2301) << "Failed sparsification with 0123->3210 nesting order change and 0123->2301 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3210>(dense_segs, lam2310);
  std::string ref_output_3210_2310 = stream.str();
  auto sparsified_3210_2310 = sparsify_iteration_space<POL3210>(dense_segs, A, perm2310);
  stream.clear();
  stream.str("");
  kernel<POL3210>(sparsified_3210_2310, lam_prt);
  std::string test_output_3210_2310 = stream.str();
  EXPECT_EQ(test_output_3210_2310,ref_output_3210_2310) << "Failed sparsification with 0123->3210 nesting order change and 0123->2310 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";

  A.impl.set_permutation(perm0123);
  stream.clear();
  stream.str(std::string());
  kernel<POL3210>(dense_segs, lam3102);
  std::string ref_output_3210_3102 = stream.str();
  auto sparsified_3210_3102 = sparsify_iteration_space<POL3210>(dense_segs, A, perm3102);
  stream.clear();
  stream.str("");
  kernel<POL3210>(sparsified_3210_3102, lam_prt);
  std::string test_output_3210_3102 = stream.str();
  EXPECT_EQ(test_output_3210_3102,ref_output_3210_3102) << "Failed sparsification with 0123->3210 nesting order change and 0123->3102 argument order change\n" << "0123 0123 reference output: " << ref_output_0123_0123 << "\n";


  
}

TYPED_TEST_SUITE_P(SparseSparsifyTest);
template <typename T>
class SparseSparsifyTest : public ::testing::Test
{
};

TYPED_TEST_P(SparseSparsifyTest, Sparsify)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  SparseTestSparsify4<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SparseSparsifyTest, Sparsify);
