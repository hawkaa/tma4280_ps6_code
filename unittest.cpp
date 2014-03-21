/* global includes */
#include <stdio.h>

/* local includes */
extern "C" {
	/* as this is a c++ program, we need to wrap it with "C" */
	#include "ps6_common_library.h"
}

#include "gtest/gtest.h"

/* matrix data base class */
class Matrix : public ::testing::Test {

protected:
	virtual void TearDown() {
		free(b);
		free(sizes);
		free(ownership);
	}

	Real **b;
	int *sizes;
	int *ownership;
	int num_ranks;
	int m;
};

/*
 * Example matrix 1
 *
 *	1 2 3
 *  	-----
 *  	4 5 6
 *	7 8 9
 */
class Matrix3x3 : public Matrix {
	protected:
	
	virtual void SetUp() {
		b = create_real_2d_array(3, 3);
		int counter = 0;
		for (int i = 0; i < 3; ++i) {
			for(int j = 0; j < 3; ++j) {
				b[i][j] = ++counter;
			}
		}
		m = 3;
		sizes = (int*)malloc(sizeof(int) * 2);
		sizes[0] = 1;
		sizes[1] = 2;

		ownership = (int*)malloc(sizeof(int) * 3);
		ownership[0] = 0;
		ownership[1] = 1;
		ownership[2] = 1;

		num_ranks = 2;

	}

};
/*
 * Example matrix 2
 *
 * 7x7, filled out row-wise with numbers 1-49
 */
class Matrix7x7 : public Matrix {
	protected:
	
	virtual void SetUp() {
		b = create_real_2d_array(7, 7);
		int counter = 0;
		for (int i = 0; i < 7; ++i) {
			for(int j = 0; j < 7; ++j) {
				b[i][j] = ++counter;
			}
		}
		m = 7;
		sizes = (int*)malloc(sizeof(int) * 3);
		sizes[0] = 2;
		sizes[1] = 2;
		sizes[2] = 3;

		ownership = (int*)malloc(sizeof(int) * 7);
		ownership[0] = 0;
		ownership[1] = 0;
		ownership[2] = 1;
		ownership[3] = 1;
		ownership[4] = 2;
		ownership[5] = 2;
		ownership[6] = 2;

		num_ranks = 3;
	}


};

/*
 * Test ownership
 */
TEST_F(Matrix3x3, create_ownership)
{
	int *ownership_output = create_ownership(m, sizes, num_ranks);
	for (int i = 0; i < m; ++i) {
		ASSERT_EQ(ownership[i], ownership_output[i]);
	}
	free(ownership_output);
}

/*
 * Test ownership
 */
TEST_F(Matrix7x7, create_ownership)
{
	int *ownership_output = create_ownership(m, sizes, num_ranks);
	for (int i = 0; i < m; ++i) {
		ASSERT_EQ(ownership[i], ownership_output[i]);
	}
	free(ownership_output);
}

/*
 * Test serial transpose
 */
TEST(transpose, HardCoded)
{
	Real **b = create_real_2d_array(2, 2);
	Real **bt = create_real_2d_array(2, 2);
	b[0][0] = 1.0;
	b[0][1] = 2.0;
	b[1][0] = 3.0;
	b[1][1] = 4.0;
	
	transpose(bt, b, 2);
	ASSERT_FLOAT_EQ(bt[0][0], 1.0);
	ASSERT_FLOAT_EQ(bt[0][1], 3.0);
	ASSERT_FLOAT_EQ(bt[1][0], 2.0);
	ASSERT_FLOAT_EQ(bt[1][1], 4.0);

	free(b);
	free(bt);
}

TEST(transpose, Looped)
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
	transpose(bt, b, matrix_size);

	/* assert if rank 0 */
	for (i = 0; i < matrix_size; ++i) {
		for  (j = 0; j < matrix_size; ++j) {
			ASSERT_FLOAT_EQ(bt[j][i], b[i][j]);
		}
	
	}
	free(b);
	free(bt);

}

/*
 * Test create_send_buffer
 */
TEST(create_send_buffer, test1)
{
	Real **b = create_real_2d_array(2, 7);	
	
	int value = 1;
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 7; j++){
			b[i][j] = value;
			++value;
		}
	}

	/* check if correctly filled out */
	ASSERT_FLOAT_EQ(3.0, b[0][2]);
	ASSERT_FLOAT_EQ(8.0, b[1][0]);

	/* generate proper test data */
	int sizes[3] = {2,2,3};
	int *s_count = create_s_count(0, 3, sizes);
	int *s_displ = create_s_displ(0, 3, sizes);
	
	Real *send_buf = create_send_buffer(b, 7, sizes, 0, 3, s_displ, s_count);

	ASSERT_FLOAT_EQ(1, send_buf[0]);
	ASSERT_FLOAT_EQ(2, send_buf[1]);
	ASSERT_FLOAT_EQ(8, send_buf[2]);
	ASSERT_FLOAT_EQ(9, send_buf[3]);
	ASSERT_FLOAT_EQ(3, send_buf[4]);
	ASSERT_FLOAT_EQ(4, send_buf[5]);
	ASSERT_FLOAT_EQ(10, send_buf[6]);
	ASSERT_FLOAT_EQ(11, send_buf[7]);
	ASSERT_FLOAT_EQ(5, send_buf[8]);
	ASSERT_FLOAT_EQ(6, send_buf[9]);
	ASSERT_FLOAT_EQ(7, send_buf[10]);
	ASSERT_FLOAT_EQ(12, send_buf[11]);
	ASSERT_FLOAT_EQ(13, send_buf[12]);
	ASSERT_FLOAT_EQ(14, send_buf[13]);

	free(b);
	free(s_count);
	free(s_displ);
	free(send_buf);
	
	
}

/*
 * Test create_sizes
 */
TEST(create_sizes, Even)
{
	int *s = create_sizes(8, 4);

	ASSERT_EQ(2, s[0]);
	ASSERT_EQ(2, s[1]);
	ASSERT_EQ(2, s[2]);
	ASSERT_EQ(2, s[3]);

	free(s);
}

/*
 * Test create_sizes
 */
TEST(create_sizes, LessRowsThanProc)
{
	int *s = create_sizes(3,5);

	ASSERT_EQ(0, s[0]);
	ASSERT_EQ(0, s[1]);
	ASSERT_EQ(1, s[2]);
	ASSERT_EQ(1, s[3]);
	ASSERT_EQ(1, s[4]);

	free(s);
}

/*
 * Test create_sizes
 */
TEST(create_sizes, OneProc)
{
	int *s = create_sizes(4,1);

	ASSERT_EQ(4, s[0]);

	free(s);
}

/*
 * test create_sizes
 */
TEST(create_sizes, UnEven)
{
	int *s = create_sizes(7, 4);

	ASSERT_EQ(1, s[0]);
	ASSERT_EQ(2, s[1]);
	ASSERT_EQ(2, s[2]);
	ASSERT_EQ(2, s[3]);

	free(s);
}

/*
 * test create_s_count
 * Testing on a 3x3 matrix for 2 proccessors
 */
TEST(create_s_count, 3x3_2p)
{
	int sizes[2] = {1, 2};
	
	/* s_count for rank 0 */

	int *s_count_r0 = create_s_count(0, 2, sizes);
	ASSERT_EQ(1, s_count_r0[0]);
	ASSERT_EQ(2, s_count_r0[1]);
	free(s_count_r0);

	int *s_count_r1 = create_s_count(1, 2, sizes);
	ASSERT_EQ(2, s_count_r1[0]);
	ASSERT_EQ(4, s_count_r1[1]);
	free(s_count_r1);
}

/*
 * test create_s_count
 * testing on a 7x7 matrix for 3 processors
 */
TEST(create_s_count, 7x7_3p)
{
	int sizes[3] = {2, 2, 3};

	int *s_count_r0 = create_s_count(0, 3, sizes);
	ASSERT_EQ(4, s_count_r0[0]);
	ASSERT_EQ(4, s_count_r0[1]);
	ASSERT_EQ(6, s_count_r0[2]);
	free(s_count_r0);

	int *s_count_r1 = create_s_count(1, 3, sizes);
	ASSERT_EQ(4, s_count_r1[0]);
	ASSERT_EQ(4, s_count_r1[1]);
	ASSERT_EQ(6, s_count_r1[2]);
	free(s_count_r1);

	int *s_count_r2 = create_s_count(2, 3, sizes);
	ASSERT_EQ(6, s_count_r2[0]);
	ASSERT_EQ(6, s_count_r2[1]);
	ASSERT_EQ(9, s_count_r2[2]);
	free(s_count_r2);
}

TEST(create_s_displ, 3x3_2p)
{
	int sizes[2] = {1, 2};

	int *s_displ_r0 = create_s_displ(0, 2, sizes);
	ASSERT_EQ(0, s_displ_r0[0]);
	ASSERT_EQ(1, s_displ_r0[1]);
	free(s_displ_r0);

	int *s_displ_r1 = create_s_displ(1, 2, sizes);
	ASSERT_EQ(0, s_displ_r1[0]);
	ASSERT_EQ(2, s_displ_r1[1]);
	free(s_displ_r1);
}

TEST(create_s_displ, 7x7_3p)
{
	int sizes[3] = {2, 2, 3};
	
	/* rank 0 */
	int *s_displ_r0 = create_s_displ(0, 3, sizes);
	ASSERT_EQ(0, s_displ_r0[0]);
	ASSERT_EQ(4, s_displ_r0[1]);
	ASSERT_EQ(8, s_displ_r0[2]);
	free(s_displ_r0);

	/* rank 1 */
	int *s_displ_r1 = create_s_displ(1, 3, sizes);
	ASSERT_EQ(0, s_displ_r1[0]);
	ASSERT_EQ(4, s_displ_r1[1]);
	ASSERT_EQ(8, s_displ_r1[2]);
	free(s_displ_r1);

	/* rank 2 */
	int *s_displ_r2 = create_s_displ(2, 3, sizes);
	ASSERT_EQ(0, s_displ_r2[0]);
	ASSERT_EQ(6, s_displ_r2[1]);
	ASSERT_EQ(12, s_displ_r2[2]);
	free(s_displ_r2);
}

TEST_F(Matrix3x3, create_matrix_rows)
{
	

	Real **b0 = create_matrix_rows(b, 3, 0, sizes);
	ASSERT_FLOAT_EQ(1, b0[0][0]);
	ASSERT_FLOAT_EQ(2, b0[0][1]);
	ASSERT_FLOAT_EQ(3, b0[0][2]);
	free(b0);

	Real **b1 = create_matrix_rows(b, 3, 1, sizes);
	ASSERT_FLOAT_EQ(4, b1[0][0]);
	ASSERT_FLOAT_EQ(5, b1[0][1]);
	ASSERT_FLOAT_EQ(6, b1[0][2]);
	ASSERT_FLOAT_EQ(7, b1[1][0]);
	ASSERT_FLOAT_EQ(8, b1[1][1]);
	ASSERT_FLOAT_EQ(9, b1[1][2]);
	free(b1);


}

TEST(get_offset, p2)
{
	int sizes[2] = {1, 2};

	ASSERT_EQ(0, get_offset(0, sizes));
	ASSERT_EQ(1, get_offset(1, sizes));
}

TEST(get_offset, Even)
{
	int sizes[4] = {2, 2, 2, 2};

	ASSERT_EQ(0, get_offset(0, sizes));
	ASSERT_EQ(2, get_offset(1, sizes));
	ASSERT_EQ(4, get_offset(2, sizes));
	ASSERT_EQ(6, get_offset(3, sizes));
}

TEST(get_offset, p3)
{
	int sizes[3] = {3, 4, 4};

	ASSERT_EQ(0, get_offset(0, sizes));
	ASSERT_EQ(3, get_offset(1, sizes));
	ASSERT_EQ(7, get_offset(2, sizes));
}

TEST(get_offest, combined)
{
	int *sizes = create_sizes(8, 4);

	ASSERT_EQ(0, get_offset(0, sizes));
	ASSERT_EQ(2, get_offset(1, sizes));
	ASSERT_EQ(4, get_offset(2, sizes));
	ASSERT_EQ(6, get_offset(3, sizes));

	free(sizes);

}

/*
 * Main function
 * Runs all tests
 */
int
main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
  	return RUN_ALL_TESTS();
}
