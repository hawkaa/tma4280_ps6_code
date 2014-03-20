/* global includes */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "mpi.h"

/*local includes */
#include "ps6_common_library.h"


static Real
get_umax(Real **b, int n, function2D u)
{
	int i, j;
	Real umax, sum;
	Real x, y;
	umax = 0.0;
	for (i = 1; i < n; ++i) {
		for (j = 1; j < n; ++j) {
			x = (Real)(j) / (Real)(n);
			y = (Real)(i) / (Real)(n);
			sum = fabs((*u)(x, y) - b[i-1][j-1]);
			if (sum > umax) {
				umax = sum;
			}
		}
	}
	return umax;

}

Real*
get_diagonal(m, n)
{
	int i;
	Real pi;
	Real *d;
	
	pi = 4.0 * atan(1.0);

	d = createRealArray(m);
	for (i = 0; i < m; ++i) {
		d[i] = 2.0 * (1.0 - cos((i + 1) * pi / (Real)n));
	}
	return d;
}


/*
 * Poisson solver.
 * Will return max error
 * Should be called from all processes
 * Only rank 0 will have valid result
 */
Real
poisson(int n, function2D f, function2D u)
{
	
	/* vector and matrix structures */
	Real *diagonal;
	Real **b;
	Real **bt;
	Real *z;
	

	Real pi, h, umax;
	int i, j,  m, nn;

	
	Real x, y;
	
	/* the total number of grid points in each spatial direction is (n+1) */
	/* the total number of degrees-of-freedom in each spatial direction is (n-1) */
	/* this version requires n to be a power of 2 */
	
	m = n - 1;
	nn = 4 * n;
	
	b = createReal2DArray (m,m);
	bt = createReal2DArray (m,m);
	z = createRealArray (nn);
	

	h = 1.0 / (Real)n;
	
	diagonal = get_diagonal(m, n);

	for (i = 0; i < m; ++i) {
		for (j = 0; j < m; ++j) {
			x = (Real)(j+1) / (Real)(n);
			y = (Real)(i+1) / (Real)(n);
			b[i][j] = h * h * (*f)(x, y);
		}
	}

	for (i = 0; i < m; ++i) {
		fst_(b[i], &n, z, &nn);
	}
	
	transpose(bt, b, m);
	
	for (i = 0; i < m; ++i) {
		fstinv_(bt[i], &n, z, &nn);
	}

	/* step 2 */
	for (i = 0; i < m; ++i) {
		for (j = 0; j < m; ++j) {
			bt[i][j] = bt[i][j]/(diagonal[i]+diagonal[j]);
		}
	}

	/* step 3 */
	for (i = 0; i < m; ++i) {
		fst_(bt[i], &n, z, &nn);
	}
	
	transpose (b, bt, m);
	
	for (i = 0; i < m; ++i) {
		fstinv_(b[i], &n, z, &nn);
	}
	
	return get_umax(b, n, u);
}





static void
print2dArray(Real** arr, int rows, int cols)
{
	int i, j;
	for(i = 0; i <rows; i++){
		for(j = 0; j < cols; j++){
			printf("%f   ", arr[i][j]);
		}
		printf("\n");
	}
}

static void
printArr(Real* arr, int size)
{
	int i;
	for(i = 0; i < size; i++){
		printf("%f   ", arr[i]);
	}
	printf("\n");
}



/*
 * Transpose function
 * Only rank 0 will have a valid result
 */
void
transpose(Real **bt, Real **b, int m)
{
	#ifdef HAVE_MPI
	/* spre ut matrise på alle prosesser */
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int* sizes = create_SIZES(m, size);
	Real** b_partial = createReal2DArray(sizes[rank], m);
	if(rank == 0){	
		int i;
		for(i = 1; i < size; i++){
			b_partial = get_matrix_rows(b, m, i, sizes);
			MPI_Send(&(b_partial[0][0]), m*sizes[i], MPI_DOUBLE, i , 100, MPI_COMM_WORLD);
		}
		b_partial = get_matrix_rows(b, m , 0, sizes);
	} else{
		MPI_Recv(&(b_partial[0][0]), m*sizes[rank], MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}	

	/* kall til funksjon som antar at matrisen allerede er spredt */
	int* s_count = create_Scount(rank, size, sizes);
	int* s_displ = create_Sdispl(rank, size, sizes);

	Real* send_buf = create_Send_buf(b_partial, rank, size, sizes, m, s_displ, s_count);
	Real* recv_buf = (Real*)malloc(sizeof(Real)*m*sizes[rank]);
	
	MPI_Alltoallv(send_buf, s_count, s_displ, MPI_DOUBLE, recv_buf, s_count, s_displ, MPI_DOUBLE, MPI_COMM_WORLD);
	
	Real** partial_trans = create_partial_transposed(recv_buf, m, rank, sizes);

	//if(rank == 0) print2dArray(partial_trans, sizes[rank], m);		

	/* samle sammen på p0 og returner */
  	if(rank == 0){	
		int i;
		for(i = 1; i < size; i++){
			MPI_Recv(&(bt[get_offset(i, sizes)][0]), m*sizes[i], MPI_DOUBLE, i, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		memcpy(bt[0], partial_trans[0], sizeof(Real) * m * sizes[0]);
	} else{
		MPI_Send(&(partial_trans[0][0]), m*sizes[rank], MPI_DOUBLE, 0 , 100, MPI_COMM_WORLD);
	}	

	#else
	int i, j;
  	for (j=0; j < m; j++) {
    		for (i=0; i < m; i++) {
      			bt[j][i] = b[i][j];
    		}
  	}
	#endif
}

Real *createRealArray (int n)
{
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

Real
**createReal2DArray(int n1, int n2)
{
	int i, n;
	Real **a;
	a    = (Real **)malloc(n1   *sizeof(Real *));
	a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
	for (i=1; i < n1; i++) {
	  a[i] = a[i-1] + n2;
	}
	n = n1*n2;
	memset(a[0],0,n*sizeof(Real));
	return (a);
}

/*  */
int
*create_SIZES(int num_rows, int num_ranks)
{
	int* sizes_arr = (int*)malloc(sizeof(int)*num_ranks);
	int i;
	/* find number of rows per process */
	int num_p_proc = num_rows/num_ranks;
	/* possible rest rows */
	int num_p_proc_r = num_rows%num_ranks;	
	/* node 0 gets less work */
	sizes_arr[0] = num_p_proc;
	for(i = num_ranks-1; i> 0; i--){
		if(num_rows <= 0){
			/* if we do not have any rows left, give zero to the rest */
			sizes_arr[i] = 0;
		} else if(num_p_proc_r == 0){
			/* if no more rows remains, give the right number of rows to process i */
			sizes_arr[i] = num_p_proc;
		} else{
			/* if we have remaining rows, give these to the last nodes */
			sizes_arr[i] = num_p_proc + 1;
			num_p_proc_r--; 
		}
		num_rows--;	
	}
	/* return sizes array */	
	return sizes_arr;
		
}

/* called by each process to know the number of elements to send to other processes */
int
*create_Scount(int current_rank, int num_ranks, int* sizes)
{
	int i;
	int* s_count = (int*)malloc(sizeof(int)*num_ranks);
	for(i = 0; i < num_ranks; i++){
		s_count[i] = sizes[current_rank]*sizes[i];	
	}
	return s_count;
}

/* called by each process to know the displacement in the send buffer to each process */
int 
*create_Sdispl(int current_rank, int num_ranks, int* sizes)
{
	int i;
	int* s_displ = (int*)malloc(sizeof(int)*num_ranks);
	s_displ[0] = 0;
	for(i = 1; i < num_ranks; i++){
		s_displ[i] = s_displ[i-1] + (sizes[current_rank]*sizes[i-1]);	
	}
	return s_displ;
}

int*
get_ownership(int num_rows, int num_ranks)
{
	int *ownership = malloc(sizeof(int) * num_rows);
	int *sizes = create_SIZES(num_rows, num_ranks);
	
	int rank, i;
	rank = 0;
	for (i = 0; i < num_rows; ++i) {
		ownership[i] = rank;
		--sizes[rank];
		if (sizes[rank] == 0)
			++rank;
	}

	return ownership;
}

Real*
create_Send_buf(Real** owned_rows, int current_rank, int num_ranks, int* sizes, int m, int* s_displ, int* s_count)
{
	int i, j;
	int index = 0;
	int num_rows = sizes[current_rank];
	int send_buffer_size = m*num_rows;
	int base;
	int rank;
	Real* send_buf = createRealArray(send_buffer_size);
	

	int* column_ownership = get_ownership(m, num_ranks);


	int offsets[num_ranks];
	for(i = 0; i < num_ranks; ++i)
		offsets[i] = 0;
	
	for(i = 0; i < num_rows; i++){
		for(j = 0; j < m; j++) {
			/* get rank for current matric element */
			rank = column_ownership[j];
			base = s_displ[rank] + (offsets[rank]++);
			send_buf[base] = owned_rows[i][j];
		}		
	}
	return send_buf;	
}

Real**
create_partial_transposed(Real* recv_buf, int m, int current_rank, int* sizes)
{
	int i, row, col;
	int current_row_num = sizes[current_rank];
	int recv_buf_length = m*current_row_num;
	Real** bt_partial = createReal2DArray(current_row_num, m);
	for(i = 0; i < recv_buf_length; i++){
		row = i % current_row_num;
		col = i / current_row_num;
		bt_partial[row][col] = recv_buf[i];
	}
	return bt_partial;
}

int
get_offset(int current_rank, int *sizes)
{
	int offset, i;
	offset = 0;
	for (i = 0; i < current_rank; ++i) {
		offset += sizes[i];
	}
	return offset;
}

Real**
get_matrix_rows(Real** b, int m, int current_rank, int *sizes)
{
	Real** b_return;
	int offset, i, j;
	
	b_return = createReal2DArray(sizes[current_rank], m);
	offset = get_offset(current_rank, sizes);
	for (i = 0; i < sizes[current_rank]; ++i) {
		memcpy(b_return[i], b[i+offset], sizeof(Real) * m);
	}

	return b_return;
}
