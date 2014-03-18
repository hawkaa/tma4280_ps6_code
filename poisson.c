/* global includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"

/* local includes */
#include "ps6_common_library.h"

Real
f(Real x, Real y)
{
	Real pi;
	pi = 4.0 * atan(1.0);
	return 5.0 * pi * pi * sin(pi * x) * sin(2.0 * pi * y);
}

Real
u(Real x, Real y)
{
	Real pi;
	pi = 4.0 * atan(1.0);
	return sin(pi * x) * sin(2.0 * pi * y);
}

void
printArr(int* arr, int size)
{
	int i;
	for(i = 0; i < size; i++){
		printf("%d   ", arr[i]);
	}
	printf("\n");
}

int
main(int argc, char** argv)
{
	int n, i, j, rank, size, tag;	
		
	/* Test */
	int counter;

	//Real **solution;
	if (argc < 2)  {
		printf("need a problem size\n");
		return 1;
	}

	n  = atoi(argv[1]);
	
	#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	Real** b;
	Real** bt;
	
	b = createReal2DArray(n, n);
	
	counter = 1;
	for (i = 0; i < n; ++i) {
		for  (j = 0; j < n; ++j) {
			b[i][j] = counter++;
		}
	}

	
	transpose(bt, b, n);

	#endif	
	
	/*int* sizes = create_SIZES(n, size);
	printf("SIZES: ");
	printArr(sizes, size);

	int* s_count = create_Scount(rank, size, sizes);
	printf("S_count: ");
	printArr(s_count, size);

	int* s_displ = create_Sdispl(rank, size, sizes);
	printf("S_displ: ");
	printArr(s_displ, size);
	
	free(sizes);
	free(s_count);
	free(s_displ);*/
	
	#ifdef HAVE_MPI
	MPI_Finalize();
	#endif

	/*
	solution = poisson(n, *f);
	//printf("Umax: %f\n", get_umax(solution, n));
	for (i = 0; i < n - 1; ++i) {
		for (j = 0; j < n - 1; ++j) {
			printf("%f\t", solution[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	Real x, y;
	Real umax, res;
	umax = 0.0;
	for (i = 1; i < n; ++i) {
		for (j = 1; j < n; ++j) {
			x = (Real)j / (Real)(n);
			y = (Real)i / (Real)(n);
			printf("x = %f\n", x);
			res = fabs(solution[i-1][j-1] - u(x, y));
			if ( res > umax) {
				umax = res;
			}
			printf("%f\t", u(x, y));
		}
		printf("\n");
	}
	printf("\nUMAX: %f\n", umax);
	*/

	return 0;

}
