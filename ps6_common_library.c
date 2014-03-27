/*
 * PROBLEM SET 6 COMMON LIBRARY
 *
 * Contains all the logic for the FST Poisson Solver for problem set 6. The
 * functions in this module is purposedly made stateless to simplify unit
 * testing.
 *
 * The parameters keep a naming convention:
 *
 * n		problem size
 * m		matrix size (problem size - 1)
 * b_part 	a part of the problem matrix, splitted in rows
 * bt_part	as b_part, but transposed
 * rank		current process/rank (in an MPI context)
 * num_ranks	number of processes/ranks (in an MPI context)
 * sizes 	array contaning how many row each rank own
 * s_count 	value that are needed by MPI to know how much data to send to each
 *		process
 * s_displ	displacement array needed by MPI to find out where to locate the
 *		data that goes to each process
 *
 */

/* global includes */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

/*local includes */
#include "ps6_common_library.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#ifdef HAVE_OPENMP
#include "omp.h"
#endif

/*
 * INTERNAL METHODS
 */

/*
 * Print 2 dimensional array helper method
 * Used while debugging
 *
 * arr - array to print
 * rows - number of rows in array
 * cols - number of columns in array
 */
static void
print_2d_array(const Real** arr, const int rows,  const int cols)
{
	int i, j;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			printf("%f   ", arr[i][j]);
		}
		printf("\n");
	}
}

/*
 * Print 1 dimensional array helper method
 *
 * arr - array to print
 * size - array length
 */
static void
print_array(const Real* arr, const int size)
{
	int i;
	for(i = 0; i < size; i++){
		printf("%f   ", arr[i]);
	}
	printf("\n");
}

/*
 * Generate the full diagonal needed for the fst poisson solver.
 */
static Real*
create_diagonal(const int m, const int n)
{
	/* loop variable */
	int i;
	
	/* pi reference */
	Real pi = 4.0 * atan(1.0);

	/* allocate array */
	Real *d = create_real_array(m);

	/* fill in data */
	for (i = 0; i < m; ++i) {
		d[i] = 2.0 * (1.0 - cos((i + 1) * pi / (Real)n));
	}

	/* return array pointer */
	return d;
}

/*
 * Allocate real array and initialize to 0.0
 *
 * n - array size
 */
Real*
create_real_array(const int n)
{
	Real *a;
	int i;
	a = (Real *)malloc(n*sizeof(Real));

	/* initialize with OMP (if available) */
	#pragma omp parallel for schedule(static) private(i)
	for (i=0; i < n; i++) {
		a[i] = 0.0;
	}
	return (a);
}

/*
 * Allocate real two dimensional array
 *
 * n1 - number of rows
 * n2 - number of columns
 */
Real**
create_real_2d_array(int n1, int n2)
{
	int i, n;
	Real **a;
	a    = (Real **)malloc(n1   *sizeof(Real *));
	a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
	/* unable to openmp parallelize */
	for (i=1; i < n1; i++) {
	  a[i] = a[i-1] + n2;
	}
	n = n1*n2;
	memset(a[0],0,n*sizeof(Real));
	return (a);
}

/*
 * De-allocates two dimensional array
 *
 * arr - the array to free
 */
void
free_real_2d_array(Real** arr)
{
	free(arr[0]);
	free(arr);
}

/*
 * Create sizes array
 * Will calculate how many rows each rank should own.
 */
int*
create_sizes(int m, const int num_ranks)
{
	/* loop variable */
	int i;

	/* allocate array */
	int* sizes_arr = (int*)malloc(sizeof(int)*num_ranks);

	/* find number of rows per process */
	int num_p_proc = m/num_ranks;

	/* possible rest rows */
	int num_p_proc_r = m%num_ranks;	

	/* node 0 gets less work */
	sizes_arr[0] = num_p_proc;

	for(i = num_ranks-1; i> 0; i--){
		if(m <= 0){
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
		m--;	
	}

	/* return sizes array */	
	return sizes_arr;
		
}

/*
 * Create s count array
 * This function will calculate how much data the current rank should send to
 * each process. Each rank element will have the current rank size multiplied
 * with own size.
 */
int*
create_s_count(const int rank, const int num_ranks, const int* sizes)
{
	/* loop variable */
	int i;

	/* allocate array */
	int* s_count = (int*)malloc(sizeof(int) * num_ranks);

	/* can be done in parallel */
	#pragma omp parallel for schedule(static) private(i)
	for (i = 0; i < num_ranks; i++) {
		s_count[i] = sizes[rank] * sizes[i];	
	}

	/* return array */
	return s_count;
}

/*
 * Create s displacement array
 * Called by each process to know the displacement in the send buffer to each
 * process
 */
int*
create_s_displ(const int rank, const int num_ranks, const int* sizes)
{
	/* loop variable */
	int i;

	/* allocate array */
	int* s_displ = (int*)malloc(sizeof(int)*num_ranks);

	/* displacement base case, rank 0 will always be 0 */
	s_displ[0] = 0;

	/* iterating over values, next value will depend on previous */
	for (i = 1; i < num_ranks; i++) {
		s_displ[i] = s_displ[i - 1] + (sizes[rank] * sizes[i - 1]);
	}

	/* return array */
	return s_displ;
}

/*
 * Create ownership array
 * Will allocate an array of m size that indicates which rank owns which row
 */
int*
create_ownership(const int m, const int* sizes, const int num_ranks)
{
	/* loop variable */
	int i;
	
	/* allocate array */
	int *ownership = malloc(sizeof(int) * m);
	
	/* start with rank 0 */
	int rank = 0;

	/* initiate row counter to size of rank 0 */
	int c = sizes[0];

	for (i = 0; i < m; ++i) {
		/* set ownership */
		ownership[i] = rank;

		/* decrease counter */
		--c;

		if (c == 0) {
			/* no more rows belongs to current rank, increment */
			++rank;
			c = sizes[rank];
		}
	}

	return ownership;
}

/*
 * Create send buffer
 * This method will allocate a send buffer array needed for the MPI Alltoallv.
 * It will place data going to the different ranks according to the s_count
 * and s_displ vector.
 */
Real*
create_send_buffer(Real **b_part, const int m, const int *sizes,
		const int rank, const int num_ranks, const int *s_displ, 
		const int *s_count)
{
	/* loop and index variables */
	int i, j;
	int base, inner_rank;

	/* number of rows to send */
	int num_rows = sizes[rank];

	/* send buffer size, number of rows times the number of columns */
	int send_buffer_size = m * num_rows;

	/* allocate send buffer */
	Real* send_buf = create_real_array(send_buffer_size);
	
	/* calculate ownership, which rank does eahc column belong to */
	int *column_ownership = create_ownership(m, sizes, num_ranks);

	/*
	 * Offset book-keeping array, incremented once a rank has received a
	 * value. Is initialized to 0
	 */
	int offsets[num_ranks];
	#pragma omp parallel for schedule(static) private(i)
	for(i = 0; i < num_ranks; ++i)
		offsets[i] = 0;
	

	/* fill the send buffer */
	for(i = 0; i < num_rows; i++){
		for(j = 0; j < m; j++) {
			/* get rank for current matrix row */
			inner_rank = column_ownership[j];

			/* calculate base address in send buffer, using
			displacement and offset. offset increased to keep track
			of that the position is taken */
			base = s_displ[inner_rank] + (offsets[inner_rank]++);

			/* save the actual value in the send buffer */
			send_buf[base] = b_part[i][j];
		}		
	}

	/* deallocate array */
	free(column_ownership);
	return send_buf;	
}

/*
 * Reconstruct partial from receive buffer.
 * This function will (in parallel) reconstruct a partial matrix from a receive
 * buffer received by MPI_Alltoallv.
 */
void
reconstruct_partial_from_receive_buffer(Real** b_part,
		const Real* receive_buffer, const int m, const int *sizes,
		const int rank)
{
	int i, row, col;
	int current_row_num = sizes[rank];
	int recv_buf_length = m * current_row_num;

	/* each element maps to a row and column based on integer division */
	#pragma omp parallel for schedule(static) private(i, row, col)
	for(i = 0; i < recv_buf_length; i++){
		row = i % current_row_num;
		col = i / current_row_num;
		b_part[row][col] = receive_buffer[i];
	}
}

/*
 * Get offset
 * Will calculate the row offset for the given rank.
 *
 * EXAMPLE
 * If we have three ranks with 10 rows each, each rank will have the following
 * offsets:
 *	r0 - 0
 *	r1 - 10
 *	r2 - 20
 */
int
get_offset(const int rank, const int *sizes)
{
	int offset, i;
	offset = 0;
	for (i = 0; i < rank; ++i) {
		offset += sizes[i];
	}
	return offset;
}

Real**
create_part_matrix(Real** b, const int m, const int rank, const int *sizes)
{
	Real** b_return;
	int offset, i, j;
	
	b_return = create_real_2d_array(sizes[rank], m);
	offset = get_offset(rank, sizes);
	for (i = 0; i < sizes[rank]; ++i) {
		memcpy(b_return[i], b[i+offset], sizeof(Real) * m);
	}

	return b_return;
}

/*
 * Transpose partial matrix 
 * This function will set up send- and receive buffers and receive a transposed
 * partial matrix 
 */
void
transpose_part(Real **bt_part, Real **b_part, int m, int *sizes, int rank, int num_ranks, int* s_displ, int* s_count)
{

	Real* send_buf = create_send_buffer(b_part, m, sizes, rank, num_ranks,
					s_displ, s_count);

	Real* recv_buf = (Real*)malloc(sizeof(Real) * m * sizes[rank]);
	
	MPI_Alltoallv(send_buf, s_count, s_displ, MPI_DOUBLE, recv_buf, 
				s_count, s_displ, MPI_DOUBLE, MPI_COMM_WORLD);

	reconstruct_partial_from_receive_buffer(bt_part, recv_buf, m,
					sizes, rank);
	free(send_buf);
	free(recv_buf);
		
}

/*
 * ENTRY METHODS
 */

/*
 * Poisson solver.
 * Will return max error
 * Should be called from all processes
 * Only rank 0 will have valid result
 */
Real
poisson_parallel(int n, function2D f, function2D u)
{
	
	/* vector and matrix structures */

	int i, j;
	Real x, y;
	int m = n - 1;
	int nn = 4 * n;

	#include <omp.h>
	/*
	 * omp_get_max_threads() - tilgjengelig
	 * omp_get_num_threads() - innsiden
	 * omp_get_thread_num() - "thread rank"
	 */


	/* mpi */
	int rank, num_ranks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

	/* reference values */
	Real h = 1.0 / (Real)n;


	/* "shared" variables (equal on all ranks) */
	int *sizes = create_sizes(m, num_ranks);
	int *s_displ = create_s_displ(rank, num_ranks, sizes);
	int *s_count = create_s_count(rank, num_ranks, sizes);

	/* "local" variables (will have different values on diffrent ranks) */
	int offset = get_offset(rank, sizes);

	/* diagonal, same diagonal generated for all */
	// Arne morten, lurt å lage hele for hver prosess?
	Real* diagonal = create_diagonal(m, n);

	/* helper structure for fst */
	// Arne morten, lurt å lage denne 2d for openmp?
	Real** z = create_real_2d_array(sizes[rank], nn);

	/* allocate needed data structures */
	Real **b_part = create_real_2d_array(sizes[rank], m);
	Real **bt_part = create_real_2d_array(sizes[rank], m);

	
	/* fill out initial data */
	#pragma omp parallel for schedule(static) private(i, j, x, y) 
	for (i = 0; i < sizes[rank]; ++i) {
		//printf("[Rank %i] Global %i, Local %i\n", rank, offset + i, i);
		for (j = 0; j < m; ++j) {
			x = (Real)(j + 1) / (Real)(n);
			y = (Real)(offset + i + 1) / (Real)(n);
			b_part[i][j] = h * h * (*f)(x, y);
		}
	}

	#pragma omp parallel for schedule(static) private(i) 
	for (i = 0; i < sizes[rank]; ++i) {
		fst_(b_part[i], &n, z[i], &nn);
	}
	
	transpose_part(bt_part, b_part, m, sizes, rank, num_ranks, s_displ, s_count);
		
	#pragma omp parallel for schedule(static) private(i) 
	for (i = 0; i < sizes[rank]; ++i) {
		fstinv_(bt_part[i], &n, z[i], &nn);
	}

	/* step 2 */
	#pragma omp parallel for schedule(static) private(i, j) 
	for (i = 0; i < sizes[rank]; ++i) {
		for (j = 0; j < m; ++j) {
			bt_part[i][j] = bt_part[i][j]/(diagonal[offset + i]+diagonal[j]);
		}
	}

	/* step 3 */
	#pragma omp parallel for schedule(static) private(i) 
	for (i = 0; i < sizes[rank]; ++i) {
		fst_(bt_part[i], &n, z[i], &nn);
	}
	
	transpose_part(b_part, bt_part, m, sizes, rank, num_ranks, s_displ, s_count);
	
	#pragma omp parallel for schedule(static) private(i) 
	for (i = 0; i < sizes[rank]; ++i) {
		fstinv_(b_part[i], &n, z[i], &nn);
	}

	/* calculate u_max */
	Real sum;
	Real u_max = 0.0;
	// arne morten
	//#pragma omp parallel for schedule(static) private(i, j, x, y)
	for (i = 0; i < sizes[rank]; ++i) {
		for (j = 0; j < m; ++j) {
			x = (Real)(j + 1) / (Real)(n);
			y = (Real)(i + offset + 1) / (Real)(n);
			sum = fabs((*u)(x, y) - b_part[i][j]);
			//#pragma omp critical
			if (sum > u_max) {
				u_max = sum;
			}
		}
	}
	Real u_max_rank_0 = -1;

	MPI_Reduce(&u_max, &u_max_rank_0, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	/* clean up */
	free_real_2d_array(b_part);
	free_real_2d_array(bt_part);
	free_real_2d_array(z);
	free(sizes);
	free(diagonal);
	free(s_count);
	free(s_displ);

	return u_max_rank_0;
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
	
	b = create_real_2d_array (m,m);
	bt = create_real_2d_array (m,m);
	z = create_real_array (nn);
	

	h = 1.0 / (Real)n;
	
	diagonal = create_diagonal(m, n);

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

	Real sum;
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

	free_real_2d_array(b);
	free_real_2d_array(bt);
	free(diagonal);
	free(z);
		
	return umax;
}



void
transpose_parallel(Real **bt, Real **b, int m)
{
	#ifdef HAVE_MPI
	/* spre ut matrise på alle prosesser */

	int rank, num_ranks;	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	
	int* sizes = create_sizes(m, num_ranks);
	int* s_displ = create_s_displ(rank, num_ranks, sizes);
	int* s_count = create_s_count(rank, num_ranks, sizes);

	/* allocate partial matrices */
	Real** b_part = create_real_2d_array(sizes[rank], m);
	Real** bt_part = create_real_2d_array(sizes[rank], m);

	/* distribute */
	if(rank == 0){	
		int i;
		for(i = 1; i < num_ranks; ++i){
			b_part = create_part_matrix(b, m, i, sizes);
			MPI_Send(&(b_part[0][0]), m*sizes[i], MPI_DOUBLE, i , 100, MPI_COMM_WORLD);
		}
		b_part = create_part_matrix(b, m , 0, sizes);
	} else{
		MPI_Recv(&(b_part[0][0]), m*sizes[rank], MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	transpose_part(bt_part, b_part, m, sizes, rank, num_ranks, s_displ, s_count);

	/* samle sammen på p0 og returner */
  	if(rank == 0){	
		int i;
		for(i = 1; i < num_ranks; i++){
			MPI_Recv(&(bt[get_offset(i, sizes)][0]), m*sizes[i], MPI_DOUBLE, i, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		memcpy(bt[0], bt_part[0], sizeof(Real) * m * sizes[0]);
	} else{
		MPI_Send(&(bt_part[0][0]), m*sizes[rank], MPI_DOUBLE, 0 , 100, MPI_COMM_WORLD);
	}
	
	/* clean up */
	free_real_2d_array(b_part);
	free_real_2d_array(bt_part);
	free(sizes);
	free(s_count);
	free(s_displ);
	
	#else
	transpose(bt, b, m);
	#endif
}

void
transpose(Real **bt, Real **b, int m)
{
	int i, j;
  	for (i = 0; i < m; ++i) {
    		for (j = 0; j < m; ++j) {
      			bt[i][j] = b[j][i];
    		}
  	}

}
