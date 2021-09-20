#include "RAJA_test-base.hpp"

TEST(LoopChain, ForallExec) {
  
  double * _a = new double[100];
  double * _b = new double[100];

  using VIEW_TYPE = RAJA::View<double, RAJA::Layout<1,int,1>>;

  VIEW_TYPE a(_a,100);
  VIEW_TYPE b(_b,100);

  for(int i = 0; i  <100; i++) {
    _a[i] = 10;
  }

  auto forall1 = RAJA::make_forall<RAJA::seq_exec>(RAJA::RangeSegment(0,100). [=](auto i){b(i) = a(i) + i;};

  forall1();

  for(int i = 0; i < 100; i++) {
    EXPECT_EQ(_b[i], 10 + i);
  }

  free(_a);
  free(_b);
}
