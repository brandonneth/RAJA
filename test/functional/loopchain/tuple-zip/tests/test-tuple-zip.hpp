#include "RAJA_test-base.hpp"

template <typename EXEC_POLICY>
void TupleZipTestImpl() {
  
//single element tuples
  using namespace RAJA;
  auto t1 = camp::make_tuple(1);
  auto t2 = camp::make_tuple(2);

  auto zipped1 = RAJA::tuple_zip(RAJA::make_tuple(t1));

  EXPECT_EQ(1, tuple_len(zipped1));

  EXPECT_EQ(1, tuple_len(camp::get<0>(zipped1)));
  EXPECT_EQ(1, camp::get<0>(camp::get<0>(zipped1)));



  auto tuple_zip_slice0 = tuple_zip_slice<0>(make_tuple(t1,t2), camp::idx_seq<0,1>{});

  EXPECT_EQ(2, tuple_len(tuple_zip_slice0));
  EXPECT_EQ(1, camp::get<0>(tuple_zip_slice0));
  EXPECT_EQ(2, camp::get<1>(tuple_zip_slice0));



  auto tozip2 = RAJA::make_tuple(t1,t2); 
  auto z2 = RAJA::tuple_zip(tozip2);


  EXPECT_EQ(1, tuple_len(z2));
  
  EXPECT_EQ(2, tuple_len(camp::get<0>(z2)));
  EXPECT_EQ(1, camp::get<0>(camp::get<0>(z2)));
  EXPECT_EQ(2, camp::get<1>(camp::get<0>(z2)));


  auto t01 = make_tuple(0,1);   
  auto t23 = make_tuple(2,3);

  auto z3 = tuple_zip(make_tuple(t01, t23));

  EXPECT_EQ(2, tuple_len(z3)) << "Number of tuples in zip of 2 element tuples\n";

  auto z3_1 = camp::get<0>(z3);
  auto z3_2 = camp::get<1>(z3);

  EXPECT_EQ(2, tuple_len(z3_1)) << "ELements in first of zip of 2 element tuples\n";
  EXPECT_EQ(2, tuple_len(z3_2)) << "ELements in first of zip of 2 element tuples\n";

  EXPECT_EQ(0, camp::get<0>(z3_1)) << "2 element zip [0][0]\n";
  EXPECT_EQ(2, camp::get<1>(z3_1)) << "2 element zip [0][1]\n";
  EXPECT_EQ(1, camp::get<0>(z3_2)) << "2 element zip [1][0]\n";
  EXPECT_EQ(3, camp::get<1>(z3_2)) << "2 element zip [1][1]\n";



  auto three_threes = tuple_repeat<3>(3);
  EXPECT_EQ(3, camp::get<0>(three_threes));
  EXPECT_EQ(3, camp::get<1>(three_threes));
  EXPECT_EQ(3, camp::get<2>(three_threes));

}

TYPED_TEST_SUITE_P(TupleZipTest);
template <typename T>
class TupleZipTest : public ::testing::Test
{
};

TYPED_TEST_P(TupleZipTest, TupleZip)
{
  using EXEC_POLICY = typename camp::at<TypeParam, camp::num<0>>::type;

  TupleZipTestImpl<EXEC_POLICY>();

}

REGISTER_TYPED_TEST_SUITE_P(TupleZipTest, TupleZip);
