#include "RAJA_test-base.hpp"
#include <string>
#include <vector>


void print_perm(auto A) {
  std::cout << "PERM: ";
  for(auto i : A.impl->get_permutation()) {
    std::cout << i << " ";
  }
  std::cout << "\n";
}

template <typename EXEC_POLICY>
void SparseDemoTestImpl() {
  using namespace RAJA;
  
  camp::idx_t N = 10;
  auto A = make_random_sparse_view<double,3,COO>(N,0.1);

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
  



}

TYPED_TEST_SUITE_P(SparseDemoTest);
template <typename T>
class SparseDemoTest : public ::testing::Test
{
};

TYPED_TEST_P(SparseDemoTest, SparseExample)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;
  SparseDemoTestImpl<EXEC_POLICY>();
}

REGISTER_TYPED_TEST_SUITE_P(SparseDemoTest, SparseExample);
