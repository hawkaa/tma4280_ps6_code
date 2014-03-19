#ifndef PS6_COMMON_LIBRARY_H_
#define PS6_COMMON_LIBRARY_H_
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Type definitios
 */
typedef double Real;

typedef Real (*function2D)(Real, Real);


/* function prototypes for external fst library */
extern void
fst_(Real *v, int *n, Real *w, int *nn);

extern void
fstinv_(Real *v, int *n, Real *w, int *nn);

/* functions defined in ps6_common_library.c */

extern Real
poisson(int problem_size, function2D f);


extern Real
*createRealArray(int n);

extern Real
**createReal2DArray(int m, int n);

extern void
transpose(Real **bt, Real **b, int m);

extern int* 
create_SIZES(int num_rows, int num_ranks);

extern int*
create_Scount(int current_rank, int num_ranks, int* sizes);

extern int*
create_Sdispl(int current_rank, int num_ranks, int* sizes);

extern Real*
create_Send_buf(Real** owned_rows, int current_rank, int num_ranks, int* sizes, int m, int* s_displ, int* s_count);

extern Real**
create_partial_transposed(Real* recv_buf, int m, int current_rank, int* sizes);

extern Real**
get_matrix_rows(Real** b, int m, int current_rank,  int *sizes);


extern int
get_offset(int current_rank, int *sizes);

extern int*
get_ownership(int num_rows, int num_ranks);

#ifdef __cplusplus
}
#endif
#endif
