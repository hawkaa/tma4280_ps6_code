#ifndef PS6_COMMON_LIBRARY_H_
#define PS6_COMMON_LIBRARY_H_

/*
 * Type definitios
 */
typedef double Real;

typedef Real (*function2D)(Real, Real);


/* function prototypes for external fst library */
void
fst_(Real *v, int *n, Real *w, int *nn);

void
fstinv_(Real *v, int *n, Real *w, int *nn);

/* functions defined in ps6_common_library.c */

Real
**poisson(int problem_size, function2D f);


Real
*createRealArray(int n);

Real
**createReal2DArray(int m, int n);

void
transpose(Real **bt, Real **b, int m);


#endif
