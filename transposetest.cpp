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


TEST(belongs_to_rank, test1)
{
	int s[3] = {2,2,3};
	ASSERT_EQ(0, belongs_to_rank(0, s, 3));
	ASSERT_EQ(0, belongs_to_rank(1, s, 3));
	ASSERT_EQ(1, belongs_to_rank(2, s, 3));
	ASSERT_EQ(1, belongs_to_rank(3, s, 3));
	ASSERT_EQ(2, belongs_to_rank(4, s, 3));
	ASSERT_EQ(2, belongs_to_rank(5, s, 3));
	ASSERT_EQ(2, belongs_to_rank(6, s, 3));
}	

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

TEST(create_Send_buf, test1)
{
	int i, j;
	Real **b;
	Real *send_buf;
	b = createReal2DArray(2, 7);
	int value = 1;
	for(i = 0; i < 2; i++){
		for(j = 0; j < 7; j++){
			b[i][j] = value;
			value++;	
		}
	}
	ASSERT_FLOAT_EQ(3.0, b[0][2]);
	ASSERT_FLOAT_EQ(8.0, b[1][0]);
	int s[3] = {2,2,3};
	
	int *sC = create_Scount(0, 3, s);
	int *sD = create_Sdispl(0, 3, s);
	
	send_buf = create_Send_buf(b, 0, 3, s, 7, sD, sC);
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

TEST(create_Scount, 3x3_2p)
{
	int s[2] = {1, 2};
	int *sC1 = create_Scount(0, 2, s);
	ASSERT_EQ(1, sC1[0]);
	ASSERT_EQ(2, sC1[1]);
	int *sC2 = create_Scount(1, 2, s);
	ASSERT_EQ(2, sC2[0]);
	ASSERT_EQ(4, sC2[1]);
}

TEST(create_Scount, 7x7_3p)
{
	int s[3] = {2, 2, 3};
	int *sC1 = create_Scount(0, 3, s);
	ASSERT_EQ(4, sC1[0]);
	ASSERT_EQ(4, sC1[1]);
	ASSERT_EQ(6, sC1[2]);
	int *sC2 = create_Scount(1, 3, s);
	ASSERT_EQ(4, sC2[0]);
	ASSERT_EQ(4, sC2[1]);
	ASSERT_EQ(6, sC2[2]);
	int *sC3 = create_Scount(2, 3, s);
	ASSERT_EQ(6, sC3[0]);
	ASSERT_EQ(6, sC3[1]);
	ASSERT_EQ(9, sC3[2]);
}

TEST(create_Sdispl, 3x3_2p)
{
	int s[2] = {1, 2};
	int *sD1 = create_Sdispl(0, 2, s);
	ASSERT_EQ(0, sD1[0]);
	ASSERT_EQ(1, sD1[1]);
	int *sD2 = create_Sdispl(1, 2, s);
	ASSERT_EQ(0, sD2[0]);
	ASSERT_EQ(2, sD2[1]);
}

TEST(create_Sdispl, 7x7_3p)
{
	int s[3] = {2, 2, 3};
	
	/* proc 0 */
	int *sD0 = create_Sdispl(0, 3, s);
	ASSERT_EQ(0, sD0[0]);
	ASSERT_EQ(4, sD0[1]);
	ASSERT_EQ(8, sD0[2]);

	/* proc 1 */
	int *sD1 = create_Sdispl(1, 3, s);
	ASSERT_EQ(0, sD1[0]);
	ASSERT_EQ(4, sD1[1]);
	ASSERT_EQ(8, sD1[2]);

	/* proc 2 */
	int *sD2 = create_Sdispl(2, 3, s);
	ASSERT_EQ(0, sD2[0]);
	ASSERT_EQ(6, sD2[1]);
	ASSERT_EQ(12, sD2[2]);
}

TEST(get_matrix_rows, 3x3)
{
	Real **b = createReal2DArray(3, 3);
	b[0][0] = 1;	b[0][1] = 2;	b[0][2] = 3;
	b[1][0] = 4;	b[1][1] = 5;	b[1][2] = 6;
	b[2][0] = 7;	b[2][1] = 8;	b[2][2] = 9;
	
	int s[2] = {1, 2};

	Real **b0 = get_matrix_rows(b, 3, 0, s);
	ASSERT_FLOAT_EQ(1, b0[0][0]);
	ASSERT_FLOAT_EQ(2, b0[0][1]);
	ASSERT_FLOAT_EQ(3, b0[0][2]);

	Real **b1 = get_matrix_rows(b, 3, 1, s);
	ASSERT_FLOAT_EQ(4, b1[0][0]);
	ASSERT_FLOAT_EQ(5, b1[0][1]);
	ASSERT_FLOAT_EQ(6, b1[0][2]);
	ASSERT_FLOAT_EQ(7, b1[1][0]);
	ASSERT_FLOAT_EQ(8, b1[1][1]);
	ASSERT_FLOAT_EQ(9, b1[1][2]);


}

TEST(get_offset, p2)
{
	int s[2] = {1, 2};
	ASSERT_EQ(0, get_offset(0, s));
	ASSERT_EQ(1, get_offset(1, s));
}

TEST(get_offset, p3)
{
	int s[3] = {3, 4, 4};
	ASSERT_EQ(0, get_offset(0, s));
	ASSERT_EQ(3, get_offset(1, s));
	ASSERT_EQ(7, get_offset(2, s));
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
