/* global includes */
#include <stdio.h>
#ifdef HAVE_MPI
#include <mpi.h>
#endif
/* local includes */
extern "C" {
#include "ps6_common_library.h"
}
#include "gtest/gtest.h"


TEST(TransposeSingleProc, HardCoded)
{
	Real **b, **bt;
	b = createReal2DArray(2, 2);
	bt = createReal2DArray(2, 2);
	b[0][0] = 1.0;
	b[0][1] = 2.0;
	b[1][0] = 3.0;
	b[1][1] = 4.0;
	
	transpose(bt, b, 2);

	ASSERT_FLOAT_EQ(bt[0][0], 1.0);
	ASSERT_FLOAT_EQ(bt[0][1], 3.0);
	ASSERT_FLOAT_EQ(bt[1][0], 2.0);
	ASSERT_FLOAT_EQ(bt[1][1], 4.0);
}

TEST(TransposeSingleProc, Looped)
{
	
	int matrix_size;

	int counter, i, j;
	matrix_size = 100;
	Real **b, **bt;
	b = createReal2DArray(matrix_size, matrix_size);
	bt = createReal2DArray(matrix_size, matrix_size);

	counter = 1;
	for (i = 0; i < matrix_size; ++i) {
		for  (j = 0; j < matrix_size; ++j) {
			b[i][j] = counter++;
		}
	}

	transpose(bt, b, matrix_size);

	for (i = 0; i < matrix_size; ++i) {
		for  (j = 0; j < matrix_size; ++j) {
			ASSERT_FLOAT_EQ(bt[j][i], b[i][j]);
		}

	}

}

TEST(create_SIZES, Even)
{
	int *s = create_SIZES(8, 4);
	ASSERT_EQ(2, s[0]);
	ASSERT_EQ(2, s[1]);
	ASSERT_EQ(2, s[2]);
	ASSERT_EQ(2, s[3]);
}

TEST(create_SIZES, LessRowsThanProc)
{
	int *s = create_SIZES(3,5);
	ASSERT_EQ(0, s[0]);
	ASSERT_EQ(0, s[1]);
	ASSERT_EQ(1, s[2]);
	ASSERT_EQ(1, s[3]);
	ASSERT_EQ(1, s[4]);
}

TEST(create_SIZES, OneProc)
{
	int *s = create_SIZES(4,1);
	ASSERT_EQ(4, s[0]);
}

TEST(create_SIZES, UnEven)
{
	int *s = create_SIZES(7, 4);
	ASSERT_EQ(1, s[0]);
	ASSERT_EQ(2, s[1]);
	ASSERT_EQ(2, s[2]);
	ASSERT_EQ(2, s[3]);
}

int
main(int argc, char** argv)
{
	int ret_val, rank;
	MPI::Init();
	::testing::InitGoogleTest(&argc, argv);

	#ifdef HAVE_MPI
	rank = MPI::COMM_WORLD.Get_rank();	
	#else
	rank = 0;
	#endif
	if (rank == 0) {
		ret_val = RUN_ALL_TESTS();
	}


	MPI::Finalize();
	return ret_val;
}
