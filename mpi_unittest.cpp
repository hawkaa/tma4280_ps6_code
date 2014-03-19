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

static int rank;


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
	if (rank != 0)
		return;
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

	if (rank != 0)
		return;

	for (i = 0; i < matrix_size; ++i) {
		for  (j = 0; j < matrix_size; ++j) {
			ASSERT_FLOAT_EQ(bt[j][i], b[i][j]);
		}

	}

}
int
main(int argc, char** argv)
{
	int ret_val;

	MPI::Init();
	::testing::InitGoogleTest(&argc, argv);

	#ifdef HAVE_MPI
	rank = MPI::COMM_WORLD.Get_rank();	
	#else
	rank = 0;
	#endif
	
  	ret_val = RUN_ALL_TESTS();

	MPI::Finalize();
	return ret_val;
}
