#include "RAJA_test-base.hpp"
#include <string>
template <typename EXEC_POLICY>
void ApplyFuseTestImpl() {
  using namespace RAJA;
  double * _a = new double[100];
  double * _b = new double[100];
  double * _c = new double[100];

  using View1 = View<double, RAJA::Layout<1, int, 1>>;

  View1 a1(_a, 100);
  View1 b1(_b, 100);
  View1 c1(_c, 100);
  using KPol1 = KernelPolicy<
    statement::For<0,RAJA::seq_exec,
      statement::Lambda<0>
    >
  >;

  for(int i = 0; i < 100; i++) {
    _a[i] = 10;
    _b[i] = 0;
  }
  auto knl1 = RAJA::make_kernel<KPol1>(
                RAJA::make_tuple(RAJA::RangeSegment(10,20)),
                [&](auto i) {b1(i) = a1(i);}
              );
  auto knl2 = RAJA::make_kernel<KPol1>(
                RAJA::make_tuple(RAJA::RangeSegment(10,20)),
                [=](auto i) {a1(i+1) = 1;}
              );

  auto fusedKnls = fuse(knl1, knl2);

  fusedKnls();


  int c1err1 = 0;
  for(int i = 10; i < 11; i++) {
    EXPECT_EQ(_b[i], 10);
  }
  for(int i = 11; i < 20; i++) {
    EXPECT_EQ(_b[i], 1);
  }



  
  free(_a);
  free(_b);
  free(_c);



}

TYPED_TEST_SUITE_P(ApplyFuseTest);
template <typename T>
class ApplyFuseTest : public ::testing::Test
{
};

TYPED_TEST_P(ApplyFuseTest, ApplyFuse)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  ApplyFuseTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(ApplyFuseTest, ApplyFuse);
