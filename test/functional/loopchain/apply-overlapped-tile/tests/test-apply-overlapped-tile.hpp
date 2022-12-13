#include "RAJA_test-base.hpp"
#include <string>
template <typename EXEC_POLICY>
void ApplyOverlappedTileTestImpl() {
 

  double * _a = new double[16*16];
  double * _b = new double[16*16];
  double * _c = new double[16*16];
  double * _d = new double[16*16];
  
  using namespace RAJA;

  using VIEW_TYPE2 = View<double, RAJA::Layout<2>>;
  VIEW_TYPE2 a(_a,16,16);
  VIEW_TYPE2 b(_b,16,16);
  VIEW_TYPE2 c(_c,16,16);
  VIEW_TYPE2 d(_d,16,16);

  auto initialize = [=]() {
    for(int i = 0; i < 16; i++) {
      for(int j = 0; j < 16; j++) {
        a(i,j) = 1;
        b(i,j) = 10;
        c(i,j) = 20;
        d(i,j) = 30;
      }
    }
  };

  auto lambda1 = [=](auto i, auto j) {
    b(i,j) = (a(i-1,j) + a(i,j-1) + a(i,j) + a(i+1,j) + a(i,j+1));
  };
  auto lambda2 = [=](auto i, auto j) {
   c(i,j) = (b(i-1,j) + b(i,j-1) + b(i,j) + b(i+1,j) + b(i,j+1));
  };
  auto lambda3 = [=](auto i, auto j) {
    d(i,j) = (c(i-1,j) + c(i,j-1) + c(i,j) + c(i+1,j) + c(i,j+1));
  };

   auto checksum = [=]() {
    int sum = 0;
    for(int i = 0; i < 16; i++) {
      for(int j = 0; j < 16; j++) {
        sum += a(i,j);
        sum += b(i,j);
        sum += c(i,j);
        sum += d(i,j);
      }
    }
    return sum;
  };

  auto rangeSegment = make_tuple(RangeSegment(1,15), RangeSegment(1,15));

  using KPOL =  RAJA::KernelPolicy<
    RAJA::statement::For<0,RAJA::omp_parallel_for_exec,
      RAJA::statement::For<1,RAJA::seq_exec,
        RAJA::statement::Lambda<0>
      >
    >
  >;

  auto knl1 = RAJA::make_kernel<KPOL>(rangeSegment, lambda1);
  auto knl2 = RAJA::make_kernel<KPOL>(rangeSegment, lambda2);
  auto knl3 = RAJA::make_kernel<KPOL>(rangeSegment, lambda3);

  initialize();
  knl1();
  knl2();
  knl3();
  auto correct_val = checksum();

  auto otiled = overlapped_tile_no_fuse<>(knl1, knl2, knl3);

  initialize();
  otiled();

  int otiledCheck = checksum();



  auto print_array = [](auto a) {
    for(int i = 0; i < 16; i++) {
      for(int j = 0; j < 16; j++) {
        std::cout << a(i,j) << " ";
      }
      std::cout << "\n";
    }

  };
  if(correct_val != otiledCheck) {
    std::cout << "Error with overlapped tiling of 3 jacobis.\n";
    std::cout << "a:\n";
    print_array(a);
    std::cout << "b:\n";
    print_array(b);
    std::cout << "c:\n";
    print_array(c);
    std::cout << "d:\n";
    print_array(d);
  } 

  
  EXPECT_EQ(correct_val, otiledCheck);
}

TYPED_TEST_SUITE_P(ApplyOverlappedTileTest);
template <typename T>
class ApplyOverlappedTileTest : public ::testing::Test
{
};

TYPED_TEST_P(ApplyOverlappedTileTest, ApplyOverlappedTile)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  ApplyOverlappedTileTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(ApplyOverlappedTileTest, ApplyOverlappedTile);
