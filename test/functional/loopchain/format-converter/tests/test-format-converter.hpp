#include "RAJA_test-base.hpp"

template <typename EXEC_POLICY>
void FormatConverterTestImpl() {

  using namespace RAJA;
  RAJA::idx_t N = 4;
  double * _a = new double[N*N];

  using KPOL = KernelPolicy<
    statement::For<0, loop_exec,
      statement::For<1, loop_exec,
          statement::Lambda<0>
        
      >
    >
  >;
  
  using VIEW2 = View<double, RAJA::Layout<2>>;
  std::array<idx_t, 2> sizes{{N,N}}; 
  auto layout_01 = make_permuted_layout(sizes, {{0,1}}); 
  auto layout_10 = make_permuted_layout(sizes, {{1,0}}); 
  
 
  VIEW2 a(_a, layout_01); 

  auto converter_10 = permute_view(a, layout_10);
  auto converter_01 = permute_view(a, layout_01);

  //set up data, no repeating values in the array
  for(int i = 0; i < N*N; i++) {
    _a[i] = i;
  }
  
  //make sure everything checks out with the view with original layout
  for(int i = 0; i < N; i++) {
    for(int j= 0; j < N; j++) {
      EXPECT_EQ(a(i,j), j+4*i) << "Initial values of view unexpected";
    }
  }

  //run conversion
  converter_10();

  //check that the view accesses are all still the same
  for(int i = 0; i < N; i++) {
    for(int j= 0; j < N; j++) {
      EXPECT_EQ(a(i,j), j+4*i) << "View access after conversion incorrect\n";
    }
  }

  //check the underlying values of the data are in line with the new layout
  for(int i = 0; i < N; i++) {
    for(int j= 0; j < N; j++) {
      int idx = j+4*i;
      int val = i+4*j;
      EXPECT_EQ(_a[idx], val) << "underlying data after conversion incorrect\n";
    }
  }

 }

TYPED_TEST_SUITE_P(FormatConverterTest);
template <typename T>
class FormatConverterTest : public ::testing::Test
{
};

TYPED_TEST_P(FormatConverterTest, FormatConverter)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  FormatConverterTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(FormatConverterTest, FormatConverter);
