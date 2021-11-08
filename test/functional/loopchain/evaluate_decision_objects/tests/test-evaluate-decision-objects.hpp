#include "RAJA_test-base.hpp"

template <typename EXEC_POLICY>
void EvaluateDecisionObjectsTestImpl() {

  using namespace RAJA;
  RAJA::idx_t N = 512;
  double * _a = new double[N*N];
  double * _b = new double[N*N];
  double * _c = new double[N*N];
  double * _d = new double[N*N];
  double * _e = new double[N*N];

  using KPOL = KernelPolicy<
    statement::For<0, loop_exec,
      statement::For<1, loop_exec,
        statement::For<2, loop_exec,
          statement::Lambda<0>
        >
      >
    >
  >;
  
  using VIEW2 = View<double, RAJA::Layout<2>>;
  VIEW2 a(_a, N, N);
  VIEW2 b(_b, N, N);
  VIEW2 c(_c, N, N);
  VIEW2 d(_d, N, N);
  VIEW2 e(_e, N, N);

  auto initialize = [&]() {
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      a(i,j) = 1;
      b(i,j) = 2;
      c(i,j) = 0;
      d(i,j) = 5;
      e(i,j) = 0;
    }
  }
  for(int i = 0; i < N/2; i++) {
    a(i,0) += 2;
    b(i,i) *= -1;
  }
  a(0,1) = -3;
  };
 
  auto checksum = [&]() {
    double s = 0;
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        s += i*c(i,j) + j*e(i,j);
      }
    }
    return s;
  }; 
  auto lambda1 = [&](auto i0, auto i1, auto i2) {
    c(i0,i2) += a(i0,i1) * b(i1, i2);
  };

  auto lambda2 = [&](auto i0, auto i1, auto i2) {
    e(i1,i0) += d(i0,i2) * c(i2,i0);
  };

   auto print_array = [=](auto a) {
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
        std::cout << a(i,j) << " ";
      }
      std::cout << "\n";
    }

  };

  auto bounds = make_tuple(RangeSegment(0,N), RangeSegment(0,N), RangeSegment(0,N));
  auto knl1 = make_kernel<KPOL>(bounds, lambda1);
  auto knl2 = make_kernel<KPOL>(bounds, lambda2);

  auto decisions = enumerate_possible_decisions<2>(tie(c));
  initialize();
  auto evaluated = evaluate_decisions(decisions, knl1, knl2);    

  

}

TYPED_TEST_SUITE_P(EvaluateDecisionObjectsTest);
template <typename T>
class EvaluateDecisionObjectsTest : public ::testing::Test
{
};

TYPED_TEST_P(EvaluateDecisionObjectsTest, EvaluateDecisionObjects)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  EvaluateDecisionObjectsTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(EvaluateDecisionObjectsTest, EvaluateDecisionObjects);
