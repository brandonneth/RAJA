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
  auto A = make_random_sparse_view3<double>(N);

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
  


  std::cout << "Printing SparseView nonzero coordinates...\n";
  kernel<POL012>(dense_segs, lam012);
  std::cout << stream.str() << "\n\n";
  
  std::cout << "\n\n-----Example 1: Canonical-----\n\n";

  std::cout << "Printing all points in dense iteration space that access nonzeros for the access A(i0,i1,i2)...\n";
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam012);
  std::cout << stream.str() << "\n\n";



  std::cout << "Sparsifying dense sparse for nesting order (0,1,2) and access A(i0,i1,i2)...\n";
  auto sparsified_012_012 = sparsify_iteration_space_equal_dims<POL012>(dense_segs, A, perm012);
  std::cout << "\n";

  std::cout << "Printing all points in the sparsified iteration space...\n";
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(sparsified_012_012, lam_prt);
  std::cout << stream.str() << "\n";  

  std::cout << "\n\n-----Example 2: Modified Access-----\n\n";
  std::cout << "Printing all points in dense iteration space that access nonzeros for the access A(i2,i0,i1)...\n";
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(dense_segs, lam201);
  std::cout << stream.str() << "\n\n";

  std::cout << "Sparsifying dense sparse for nesting order (0,1,2) and access A(i2,i0,i1)...\n";
  auto sparsified_012_201 = sparsify_iteration_space_equal_dims<POL012>(dense_segs, A, perm201);
  std::cout << "\n";

  std::cout << "Printing all points in the sparsified iteration space...\n";
  stream.clear();
  stream.str(std::string());
  kernel<POL012>(sparsified_012_201, lam_prt);
  std::cout << stream.str() << "\n";  


  std::cout << "\n\n-----Example 3: Modified Nesting-----\n\n";
  std::cout << "Printing all points in dense iteration space that access nonzeros in A(i0,i1,i2) with nesting order (1,0,2)...\n";
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(dense_segs, lam012);
  std::cout << stream.str() << "\n\n";

  std::cout << "Sparsifying dense space for nesting order (1,0,2) and access A(i0,i1,i2)...\n";
  auto sparsified_102_012 = sparsify_iteration_space_equal_dims<POL102>(dense_segs, A, perm012);
  std::cout << "\n";

  std::cout << "Printing all points in the sparsified iteration space...\n";
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(sparsified_102_012, lam_prt);
  std::cout << stream.str() << "\n";  

  std::cout << "\n\n-----Example 4: Modified Nesting and Access-----\n\n";
  std::cout << "Printing all points in dense iteration space that access nonzeros in A(i1,i2,i0) with nesting order (1,0,2)...\n";
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(dense_segs, lam120);
  std::cout << stream.str() << "\n\n";

  std::cout << "Sparsifying dense space for nesting order (1,0,2) and access A(i1,i2,i0)...\n";
  auto sparsified_102_120 = sparsify_iteration_space_equal_dims<POL102>(dense_segs, A, perm120);
  std::cout << "\n";

  std::cout << "Printing all points in the sparsified iteration space...\n";
  stream.clear();
  stream.str(std::string());
  kernel<POL102>(sparsified_102_120, lam_prt);
  std::cout << stream.str() << "\n";  





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
