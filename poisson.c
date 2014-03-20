/* global includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"

/* local includes */
#include "ps6_common_library.h"

/* Taken from the lecturers common.c library */
Real WallTime ()
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

Real*
Max_n_Min(Real* arr, int num_elem){
	int i;
	Real max, min;
	Real* r_arr;
	r_arr = (Real*)malloc(sizeof(Real)*2);
	max = min = arr[0];
	for(i = 1; i < num_elem; ++i){
		if(arr[i] > max) max = arr[i];
		if(arr[i] < min) min = arr[i];
	}
	r_arr[0] = max;
	r_arr[1] = min;
	return r_arr;

}

int
main(int argc, char** argv)
{
	int n, i, j, rank, size, tag, num_of_runs;	
	Real umax, t1, t2;
	Real* wtimes;	
	
	if (argc < 3)  {
		printf("need a problem size and number of runs\n");
		return 1;
	}

	n  = atoi(argv[1]);
	num_of_runs = atoi(argv[2]);
	
	#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif		
	 
	if(rank == 0) wtimes = (Real*)malloc(sizeof(Real)*num_of_runs);
	for(i = 0; i < num_of_runs; ++i){
		t1 = WallTime();
		umax = poisson(n, f ,u);	
		t2 = WallTime();
		if(rank == 0){
			printf("run %d, time=%f, umax=%.10f\n", i, t2-t1, umax);
			wtimes[i] = t2- t1;
		}
	}
	

	/* find avarage walltime */
	if(rank == 0){
		int i;
		Real sum;
		Real max, min;
		max = min = wtimes[0];
		sum = 0.0;
		for(i = 1; i < num_of_runs; ++i){
			if(wtimes[i] > max) max = wtimes[i];
			if(wtimes[i] < min) min = wtimes[i];
			sum += wtimes[i];
		}
		sum = sum - max - min;
		printf("min: %f\n", min);
		printf("max: %f\n", max);
		printf("avarage time: %f\n", sum/(num_of_runs-2));
	}
	
	#ifdef HAVE_MPI
	MPI_Finalize();
	#endif
	
	/*
	Real** b;
	Real** bt;
	*/
	/* Test */
	/*
	int counter;

	b = createReal2DArray(n, n);
	
	counter = 1;
	for (i = 0; i < n; ++i) {
		for  (j = 0; j < n; ++j) {
			b[i][j] = counter++;
		}
	}

	transpose(bt, b, n);
        */

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
