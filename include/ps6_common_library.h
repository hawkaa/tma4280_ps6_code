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
**poisson(int problem_size, function2D f);


extern Real
*createRealArray(int n);

extern Real
**createReal2DArray(int m, int n);

extern void
transpose(Real **bt, Real **b, int m);

#ifdef __cplusplus
}
#endif
#endif
