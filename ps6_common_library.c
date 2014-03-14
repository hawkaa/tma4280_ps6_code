/* global includes */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

/*local includes */
#include "ps6_common_library.h"

Real
**poisson(int problem_size, function2D f)
{
	Real *diag, **b, **bt, *z;
	Real pi, h, umax;
	int i, j, n, m, nn;
	
	/* the total number of grid points in each spatial direction is (n+1) */
	/* the total number of degrees-of-freedom in each spatial direction is (n-1) */
	/* this version requires n to be a power of 2 */
	
	n = problem_size;
	m = n - 1;
	nn = 4 * n;
	
	diag = createRealArray (m);
	b = createReal2DArray (m,m);
	bt = createReal2DArray (m,m);
	z = createRealArray (nn);
	
	h = 1.0 / (Real)n;
	pi = 4.0 * atan(1.0);
	
	for (i=0; i < m; i++) {
		diag[i] = 2.0 * (1.0 - cos((i + 1) * pi / (Real)n));
	}
	Real x, y;
	for (j=0; j < m; j++) {
		for (i=0; i < m; i++) {
			x = (Real)(j+1) / (Real)(n);
			y = (Real)(j+1) / (Real)(n);
			b[j][i] = h * h * (*f)(x, y);
		}
	}
	for (j=0; j < m; j++) {
		fst_(b[j], &n, z, &nn);
	}
	
	transpose(bt, b, m);
	
	for (i=0; i < m; i++) {
	  fstinv_(bt[i], &n, z, &nn);
	}

	/* step 2 */
	for (j=0; j < m; j++) {
	  for (i=0; i < m; i++) {
	    bt[j][i] = bt[j][i]/(diag[i]+diag[j]);
	  }
	}

	/* step 3 */
	for (i=0; i < m; i++) {
	  fst_(bt[i], &n, z, &nn);
	}
	
	transpose (b,bt,m);
	
	for (j=0; j < m; j++) {
	  fstinv_(b[j], &n, z, &nn);
	}
	
	return b;

	umax = 0.0;
	for (j=0; j < m; j++) {
	  for (i=0; i < m; i++) {
	    if (b[j][i] > umax) umax = b[j][i];
	  }
	}


	printf (" umax = %e \n",umax);
	
	return 0;
}

Real
get_umax_old(Real **solution, int problem_size)
{
	int i, j;
	Real umax;
	umax = 0.0;
	for (i = 0; i < problem_size - 1; ++i) {
		for (j = 0; j < problem_size - 1; ++j) {
			if (solution[i][j] > umax) {
				umax = solution[i][j];
			}
		}
	}
	return umax;

}



void transpose (Real **bt, Real **b, int m)
{
  int i, j;
  for (j=0; j < m; j++) {
    for (i=0; i < m; i++) {
      bt[j][i] = b[i][j];
    }
  }
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
