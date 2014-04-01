/* global includes */
#include <stdio.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/* local includes */
extern "C" {
	/* as this is a c++ program, we need to wrap it with "C" */
	#include "ps6_common_library.h"
}

#include "gtest/gtest.h"

/* save current rank */
static int rank;

/*
 * Unitest to check that the parallel matrix transpose function runs.
 * This one is a simple hard-coded 2x2 matrix.
 */
TEST(transpose_parallel, HardCoded)
{
	/* matrix buffers */
	Real **b = create_real_2d_array(2, 2);
	Real **bt = create_real_2d_array(2, 2);

	/* values */
	b[0][0] = 1.0;
	b[0][1] = 2.0;
	b[1][0] = 3.0;
	b[1][1] = 4.0;
	
	/* run transpose function */
	transpose_parallel(bt, b, 2);

	if (rank == 0) {
		/* assert if rank 0 */
		ASSERT_FLOAT_EQ(bt[0][0], 1.0);
		ASSERT_FLOAT_EQ(bt[0][1], 3.0);
		ASSERT_FLOAT_EQ(bt[1][0], 2.0);
		ASSERT_FLOAT_EQ(bt[1][1], 4.0);
	}
	free_real_2d_array(b);
	free_real_2d_array(bt);
}

/*
 * Unitest to check that the parallel matrix transpose function runs.
 * this one creates a 100x100 matrix an checks the 
 */
TEST(Transpose, Looped)
{
	
	/* loop variables */
	int i, j;
	
	/* matrix size */
	int matrix_size = 100;

	/* matrix buffers */
	Real **b = create_real_2d_array(matrix_size, matrix_size);
	Real **bt = create_real_2d_array(matrix_size, matrix_size);

	int counter = 0;
	for (i = 0; i < matrix_size; ++i) {
		for  (j = 0; j < matrix_size; ++j) {
			b[i][j] = ++counter;
		}
	}
	
	/* run the transpose function */
	transpose_parallel(bt, b, matrix_size);

	if (rank == 0) {
		/* assert if rank 0 */
		for (i = 0; i < matrix_size; ++i) {
			for  (j = 0; j < matrix_size; ++j) {
				ASSERT_FLOAT_EQ(bt[j][i], b[i][j]);
			}
	
		}
	}
	free_real_2d_array(b);
	free_real_2d_array(bt);

}

/*
 * Main function
 * Program should only run if MPI is enabled
 */
int
main(int argc, char** argv)
{
	#ifdef HAVE_MPI

	int ret_val;

	/* init mpi and gtest */
	MPI::Init();
	::testing::InitGoogleTest(&argc, argv);

	/* save rank */
	rank = MPI::COMM_WORLD.Get_rank();	
	
	/* run all tests */
  	ret_val = RUN_ALL_TESTS();

	/* finalize MPI */
	MPI::Finalize();
	
	exit(ret_val);

	#else

	printf("MPI must be enabled to run these tests.\n");
	exit(1);

	#endif
	
}
