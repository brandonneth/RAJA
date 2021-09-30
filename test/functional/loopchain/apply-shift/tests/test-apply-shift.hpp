#include "RAJA_test-base.hpp"
#include <string>
template <typename EXEC_POLICY>
void ApplyShiftTestImpl() {
  
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

  for(int i = 0; i < 100; i++) {
    _a[i] = i;
    _b[i] = 0;
  }
  auto knl1 = RAJA::make_kernel<KERNEL_POL1>(
                RAJA::make_tuple(RAJA::RangeSegment(10,20)),
                [&](auto i) {b1(i) = a1(i) + i;}
              );

  auto knl1Shifted = RAJA::shift(knl1, 5);

  knl1Shifted();

  for(int i = 0; i < 10; i++) {
    EXPECT_EQ(b1(i), 0);
   }

  for(int i = 10; i < 20; i++) {
    EXPECT_EQ(b1(i), a1(i) + i);
  }

  for(int i = 20; i < 100; i++) {
    EXPECT_EQ(b1(i), 0);
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

  using namespace RAJA;
  auto knl2 = RAJA::make_kernel<KERNEL_POL2>(
    make_tuple(RangeSegment(0,10),RangeSegment(0,10)),
    [=] (auto i, auto j) {
      b2(i,j) = a2(i,j) + 2;
    });

  auto shifted2 = RAJA::shift(knl2, 1,-1);

  for(int i = 0; i < 100; i++) {
    _a[i] = i;
    _b[i] = 0;
  }

  shifted2();


  for(int i = 0; i < 100; i++) {
    EXPECT_EQ(_b[i], i + 2);
  }



  
  free(_a);
  free(_b);



}

TYPED_TEST_SUITE_P(ApplyShiftTest);
template <typename T>
class ApplyShiftTest : public ::testing::Test
{
};

TYPED_TEST_P(ApplyShiftTest, ApplyShift)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  ApplyShiftTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(ApplyShiftTest, ApplyShift);
