#include "RAJA_test-base.hpp"

template <typename EXEC_POLICY>
void KernelExecTestImpl() {
  
  double * _a = new double[100];
  double * _b = new double[100];

  using VIEW_TYPE = RAJA::View<double, RAJA::Layout<1, int, 1>>;

  VIEW_TYPE a1(_a, 100);
  VIEW_TYPE b1(_b, 100);

  for(int i = 0; i < 100; i++) {
    _a[i] = 10;
  }

  using KERNEL_POL1 = RAJA::KernelPolicy<
    RAJA::statement::For<0,EXEC_POLICY,
      RAJA::statement::Lambda<0>
    >
  >;
  auto knl1 = RAJA::make_kernel<KERNEL_POL1>(camp::make_tuple(RAJA::RangeSegment(0,100)), [=](auto i) {b1(i) = a1(i) + i;});

  knl1();

  for(int i = 0; i < 100; i++) {
    EXPECT_EQ(_b[i], 10+i);
  }



  for(int i = 0; i < 100; i++) {
    _a[i] = i;
    _b[i] = -100;
  }

  using VIEW_TYPE2 = RAJA::View<double, RAJA::Layout<2>>;

  VIEW_TYPE2 a2(_a, 10, 10);
  VIEW_TYPE2 b2(_b, 10, 10);

  using KERNEL_POL2 = RAJA::KernelPolicy<
    RAJA::statement::For<0,RAJA::seq_exec,
      RAJA::statement::For<1,EXEC_POLICY,
        RAJA::statement::Lambda<0>
      >
    >
  >;

  auto knl2 = RAJA::make_kernel<KERNEL_POL2>(
    camp::make_tuple(
      RAJA::RangeSegment(1,9),
      RAJA::RangeSegment(0,10)),
    [=] (auto i, auto j) {
      b2(i,j) = (a2(i-1,j) + a2(i,j) + a2(i+1,j)) / 3;
    }
  );

  knl2();

  for(int i = 1; i < 9; i++) {
    for(int j = 0; j < 10; j++) {
      ASSERT_EQ(b2(i,j), (a2(i-1,j) + a2(i,j) + a2(i+1,j)) / 3);
    }
  }



  free(_a);
  free(_b);



}

TYPED_TEST_SUITE_P(KernelExecTest);
template <typename T>
class KernelExecTest : public ::testing::Test
{
};

TYPED_TEST_P(KernelExecTest, KernelExec)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  KernelExecTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(KernelExecTest, KernelExec);
