/* include guard */
#ifndef PS6_COMMON_LIBRARY_H_
#define PS6_COMMON_LIBRARY_H_

/* if c++, we need to export header as a C library */
#ifdef __cplusplus
extern "C" {
#endif

/*
 * TYPE DEFINITIONS
 */

/* floating point type used for all calculations */
typedef double Real;

/* function on the format f(x, y) */
typedef Real (*function2D)(Real, Real);


/*
 * FAST SINE TRANSFORM LIBRARY
 * Function prototypes for fst.f library
 */

/* fast sine transform */
extern void
fst_(Real *v, int *n, Real *w, int *nn);


/* fast inverse sine transform */
extern void
fstinv_(Real *v, int *n, Real *w, int *nn);


/*
 * ENTRY METHODS
 * Here comes the methods shat should be called from outside
 * Function names suffixed with "_parallel" should be called from all
 * processors in an MPI setting, but only rank 0 will have a valid result.
 * Function documentation resides in ps6_common_library.c
 */

extern Real
poisson_parallel(int problem_size, Real *time, function2D f, function2D u);

extern Real
poisson(int problem_size, function2D f, function2D u);

extern void
transpose(Real **bt, Real **b, int m);

extern void
transpose_parallel(Real **bt, Real **b, int m);

/*
 * UTILITY METHODS
 * These functions should be used from other modules, but are only thought
 * of as utility, and has no side effects regarding use of OpenMP and MPI.
 * Documentation of each function resides in ps6_common_library.c
 */

extern Real
wall_time(void);

extern Real
get_average(Real *arr, int n, int cutoff);
/*
 * INTERNAL METHODS
 * Helper methods for the ps6_common_library.c.
 * These are exposed in this header file for so the unit test framework
 * know their signatures.
 * Documentation of each function resides in ps6_common_library.c
 */

extern void
transpose_part(Real **bt_part, Real **b_part, int m, int *sizes, int rank, int num_ranks, int* s_displ, int* s_count);

extern Real*
create_real_array(const int n);

extern Real**
create_real_2d_array(int m, int n);

extern void
free_real_2d_array(Real **arr);

extern int* 
create_sizes(int m, const int num_ranks);

extern int*
create_s_count(const int rank, const int num_ranks, const int* sizes);

extern int*
create_s_displ(const int rank, const int num_ranks, const int* sizes);

extern int*
create_ownership(const int m, const int* sizes, const int num_ranks);

extern Real*
create_send_buffer(Real **b_part, const int m, const int *sizes,
		const int rank, const int num_ranks, const int *s_displ, 
		const int *s_count);

extern void
reconstruct_partial_from_receive_buffer(Real** b_part,
		const Real* receive_buffer, const int m, const int *sizes, 
		const int rank);

extern int
get_offset(const int rank, const int *sizes);

extern Real**
create_part_matrix(Real** b, const int m, const int current_rank,
		const int *sizes);



#ifdef __cplusplus
}
#endif
#endif
