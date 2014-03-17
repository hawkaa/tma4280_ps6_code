/* global includes */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

/*local includes */
#include "ps6_common_library.h"



/*
 * Poisson solver.
 * Should only be called from processor 0
 */
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
			x = (Real)(i+1) / (Real)(n);
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

}

/*
 * Transpose function
 * SHOULD ONLY BE CALLED FROM PROCESSOR 0
 */
void
transpose(Real **bt, Real **b, int m)
{
	/* spre ut matrise på alle prosesser */
		

	/* kall til funksjon som antar at matrisen allerede er spredt */


	/* samle sammen på p0 og returner */
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
	for(i = 1; i< num_ranks; i++){
		if(num_rows <= 0){
			/* if we do not have any rows left, give zero to the rest */
			sizes_arr[i] = 0;
		} else if(num_p_proc_r == 0){
			/* if no more rest rows, give the right number of rows to process i */
			sizes_arr[i] = num_p_proc;
		} else{
			/* if we have rest rows, give these to the last nodes */
			sizes_arr[i] = num_p_proc + 1;
			num_p_proc_r--; 
		}
		num_rows--;	
	}
	/* return sizes array */	
	return sizes_arr;
		
}
