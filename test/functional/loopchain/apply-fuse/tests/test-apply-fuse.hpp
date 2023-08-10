#include "RAJA_test-base.hpp"
#include <string>
template <typename EXEC_POLICY>
void ApplyFuseTestImpl() {
  using namespace RAJA;
  double * _a = new double[100];
  double * _b = new double[100];

  using View1 = View<double, RAJA::Layout<1, int, 1>>;

  View1 a1(_a, 100);
  View1 b1(_b, 100);
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


  for(int i = 10; i < 11; i++) {
    EXPECT_EQ(_b[i], 10);
  }
  for(int i = 11; i < 20; i++) {
    EXPECT_EQ(_b[i], 1);
  }


  using View2 = View<double, RAJA::Layout<2>>;

  View2 a2(_a, 10,10);
  View2 b2(_b, 10,10);

  using KPol2 = KernelPolicy<
    statement::For<0,seq_exec,
      statement::For<1,seq_exec,
        statement::Lambda<0>
      >
    >
  >;

  for(int i = 0; i < 100; i++) {
    _a[i] = 10;
    _b[i] = 0;
  }

  auto knl3 = RAJA::make_kernel<KPol2>(RAJA::make_tuple(RangeSegment(1,2), RangeSegment(0,9)),
   [=](auto i, auto j) {b2(i,j) = a2(i,j);}
  );

   auto knl4 = RAJA::make_kernel<KPol2>(RAJA::make_tuple(RangeSegment(1,2), RangeSegment(0,9)),
   [=](auto i, auto j) {a2(i,j+1) = 1;}
  );

  auto f2 = fuse(knl3, knl4);
  
  f2();

  for(int i = 10; i < 11; i++) {
    EXPECT_EQ(_b[i], 10);
  }
  for(int i = 11; i < 19; i++) {
    EXPECT_EQ(_b[i], 1);
  }
  
  free(_a);
  free(_b);



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
