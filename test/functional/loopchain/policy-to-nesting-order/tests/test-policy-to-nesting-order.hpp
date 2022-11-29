#include "RAJA_test-base.hpp"

template <typename EXEC_POLICY>
void PolicyToNestingOrderTestImpl() {
  
  

  using POL1 = RAJA::KernelPolicy<
    RAJA::statement::For<0,EXEC_POLICY,
      RAJA::statement::Lambda<0>
    >
  >;

  auto order1 = RAJA::policy_to_nesting_order<POL1>();

  EXPECT_EQ(order1[0], 0);
  
  using POL2 = RAJA::KernelPolicy<
    RAJA::statement::For<0,RAJA::seq_exec,
      RAJA::statement::For<1,EXEC_POLICY,
        RAJA::statement::Lambda<0>
      >
    >
  >;

  auto order2 = RAJA::policy_to_nesting_order<POL2>();
  EXPECT_EQ(order2[0], 0);
  EXPECT_EQ(order2[1], 1);
  using POL3 = RAJA::KernelPolicy<
    RAJA::statement::For<1,RAJA::seq_exec,
      RAJA::statement::For<0,EXEC_POLICY,
        RAJA::statement::Lambda<0>
      >
    >
  >;

  auto order3 = RAJA::policy_to_nesting_order<POL3>();
  EXPECT_EQ(order3[0], 1);
  EXPECT_EQ(order3[1], 0);

  using POL4 = RAJA::KernelPolicy<
    RAJA::statement::For<1,RAJA::seq_exec,
      RAJA::statement::For<0,EXEC_POLICY,
        RAJA::statement::For<2,EXEC_POLICY,
          RAJA::statement::Lambda<0>
        >
      >
    >
  >;

  auto order4 = RAJA::policy_to_nesting_order<POL4>();
  EXPECT_EQ(order4[0], 1);
  EXPECT_EQ(order4[1], 0);
  EXPECT_EQ(order4[2], 2);

  
}

TYPED_TEST_SUITE_P(PolicyToNestingOrderTest);
template <typename T>
class PolicyToNestingOrderTest : public ::testing::Test
{
};

TYPED_TEST_P(PolicyToNestingOrderTest, PolicyToNestingOrder)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  PolicyToNestingOrderTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(PolicyToNestingOrderTest, PolicyToNestingOrder);
