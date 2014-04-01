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
 * Helper methods for the ps6_common_library.c.
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
	#pragma omp parallel for private(i)
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
 * Will calculate how many rows each rank own.
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
 */
int
get_offset(const int rank, const int *sizes)
{
	int offset, i;
	offset = 0;
	/* accumulate sizes up to current rank */
	for (i = 0; i < rank; ++i) {
		offset += sizes[i];
	}
	return offset;
}

/*
 * Create partial matrix
 * Will construct a partial matrix based on b. It will use the processor rank
 * to calculate which rows belongs to the partial matrix.
 */
Real**
create_part_matrix(Real** b, const int m, const int rank, const int *sizes)
{
	int i, j;
	
	Real** b_return = create_real_2d_array(sizes[rank], m);
	int offset = get_offset(rank, sizes);
	for (i = 0; i < sizes[rank]; ++i) {
		/* copy the given rows into the partial matrix */
		memcpy(b_return[i], b[i+offset], sizeof(Real) * m);
	}

	return b_return;
}

/*
 * Transpose partial matrix 
 * This function will set up send buffer with the contents of b_part, allocate
 * a receive buffer and reconstruct the transposed matrix into bt_part.
 */
void
transpose_part(Real **bt_part, Real **b_part, int m, int *sizes, int rank,
		int num_ranks, int* s_displ, int* s_count)
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
 * Compare real function used for qsort int get_average
 */
static int
compare_real (const void * a, const void * b)
{
	if (*(Real*)a < *(Real*)b) return -1;
	if (*(Real*)a == *(Real*)b) return 0;
	if (*(Real*)a > *(Real*)b) return 1;
}

/*
 * UTILITY METHODS
 * These functions should be used from other modules, but are only thought
 * of as utility, and has no side effects regarding use of OpenMP and MPI.
 */


/*
 * Get average with cutoff
 * This funciton will return the average of a list with Real, removing 
 * 2*cutoff outliers.
 */
Real
get_average(Real *arr, int n, int cutoff)
{
	int i;
	
	/* sort the values */
	qsort(arr, n, sizeof(Real), &compare_real);	

	/* sum the elements, exclude outliers */
	Real sum = 0.0;
	for (i = cutoff; i < n - cutoff; ++i) {
		sum += arr[i];
	}

	/* calculate and return average */
	return sum / (n - 2 * cutoff);
}

/*
 * Wall time
 * Taken from the lecturers common.c library
 */
Real
wall_time ()
{
	#ifdef HAVE_MPI
 	return MPI_Wtime();
	#elif defined(HAVE_OPENMP)
 	return omp_get_wtime();
	#else
 	struct timeval tmpTime;
 	gettimeofday(&tmpTime,NULL);
 	return tmpTime.tv_sec + tmpTime.tv_usec/1.0e6;
	#endif
}


/*
 * ENTRY METHODS
 * Here comes the methods that should be called from outside
 * Function names suffixed with "_parallel" should be called from all
 * processors in an MPI setting, but only rank 0 will have a valid result.
 */

/*
 * Parallel poisson solver.
 * Will return max error if reference function u is present
 * Should be called from all processes
 * Only rank 0 will have valid result
 */
Real
poisson_parallel(int n, Real *time, function2D f, function2D u)
{
	
	/* vector and matrix structures */
	int i, j;
	Real x, y;
	int m = n - 1;
	int nn = 4 * n;
	Real t1, t2;



	/* mpi variables */
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
	Real* diagonal = create_diagonal(m, n);

	/* helper structure for fst */
	Real** z = create_real_2d_array(omp_get_max_threads(), nn);

	/* allocate needed data structures */
	Real **b_part = create_real_2d_array(sizes[rank], m);
	Real **bt_part = create_real_2d_array(sizes[rank], m);

	/* add a barrier to improve accuracy in timing */
	MPI_Barrier(MPI_COMM_WORLD);

	/* note start time */
	t1 = wall_time();
	
	/* fill out initial data */
	#pragma omp parallel for schedule(static) private(i, j, x, y) 
	for (i = 0; i < sizes[rank]; ++i) {
		for (j = 0; j < m; ++j) {
			/* map local indexes to x and y values */
			x = (Real)(j + 1) / (Real)(n);
			y = (Real)(offset + i + 1) / (Real)(n);
			b_part[i][j] = h * h * (*f)(x, y);
		}
		/* perform fast sine transform */
		fst_(b_part[i], &n, z[omp_get_thread_num()], &nn);
	}

	/* transpose b_part, and get a transposed partial matrix loaded into
	 bt_part */
	transpose_part(bt_part, b_part, m, sizes, rank, num_ranks, s_displ, 
			s_count);
	
	/* new openMP fork */
	#pragma omp parallel for schedule(static) private(i, j) 
	for (i = 0; i < sizes[rank]; ++i) {
		fstinv_(bt_part[i], &n, z[omp_get_thread_num()], &nn);
		for (j = 0; j < m; ++j) {
			bt_part[i][j] = bt_part[i][j]/(diagonal[offset + i]+diagonal[j]);
		}
		fst_(bt_part[i], &n, z[omp_get_thread_num()], &nn);
	}
	
	/* transpose the matrix back again, now b_part will have the interesting content */
	transpose_part(b_part, bt_part, m, sizes, rank, num_ranks, s_displ, s_count);
	
	#pragma omp parallel for schedule(static) private(i) 
	for (i = 0; i < sizes[rank]; ++i) {
		fstinv_(b_part[i], &n, z[omp_get_thread_num()], &nn);
	}
	
	
	/* register the time after all processes are done */
	MPI_Barrier(MPI_COMM_WORLD);
	t2 = wall_time();
	
	/* save time at the time pointer location (if present) */
	if (time != NULL)
		*time = t2 - t1;

	/* if we have no reference function, return invalid value */
	if (u == NULL) {
		return -1;
	}

	/* since we have a reference function, we need to calculate the maximum
	error */
	Real sum;
	Real u_max = 0.0;

	for (i = 0; i < sizes[rank]; ++i) {
		for (j = 0; j < m; ++j) {
			x = (Real)(j + 1) / (Real)(n);
			y = (Real)(i + offset + 1) / (Real)(n);
			sum = fabs((*u)(x, y) - b_part[i][j]);
			if (sum > u_max) {
				u_max = sum;
			}
		}
	}

	/* reduce umax to rank 0 */
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
 * Will return max error if reference function u is present.
 * This function is not documented very well, as it is quite similar to
 * the parallel implementation.
 */
Real
poisson(int n, Real *time, function2D f, function2D u)
{
	
	/* vector and matrix structures */
	Real *diagonal;
	Real **b;
	Real **bt;
	Real *z;
	

	Real pi, h, umax;
	int i, j,  m, nn;
	Real t1, t2;

	
	Real x, y;
	
	
	m = n - 1;
	nn = 4 * n;
	
	b = create_real_2d_array (m,m);
	bt = create_real_2d_array (m,m);
	z = create_real_array (nn);
	

	h = 1.0 / (Real)n;


	/* note start time */
	t1 = wall_time();
	
	diagonal = create_diagonal(m, n);
	Real x_left, x_right, y_top, y_bottom;
	for (i = 0; i < m; ++i) {
		for (j = 0; j < m; ++j) {
			x = (Real)(j + 1) / (Real)(n);
			y = (Real)(i + 1) / (Real)(n);
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


	t2 = wall_time();
	
	/* save time at the time pointer location (if present) */
	if (time != NULL)
		*time = t2 - t1;

	if (u == NULL) {
		return -1;
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


/*
 * Transpose
 * Will transpose b and save it to bt in paralell
 */
void
transpose_parallel(Real **bt, Real **b, int m)
{
	/* mpi */
	int rank, num_ranks;	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	
	int* sizes = create_sizes(m, num_ranks);
	int* s_displ = create_s_displ(rank, num_ranks, sizes);
	int* s_count = create_s_count(rank, num_ranks, sizes);

	/* allocate partial matrices */
	Real** b_part = create_real_2d_array(sizes[rank], m);
	Real** bt_part = create_real_2d_array(sizes[rank], m);

	/* distribute array row-wise to all ranks */
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
	
	/* transpose partial array */
	transpose_part(bt_part, b_part, m, sizes, rank, num_ranks, s_displ, s_count);

	/* gather the result back on rank 0 */
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
	
}

/*
 * Transpose
 * Saves the transposed b into bt
 */
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
