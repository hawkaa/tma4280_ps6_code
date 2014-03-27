/* global includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"

/* local includes */
#include "ps6_common_library.h"

/* constants */
#define OUTLIERS_CUTOFF 2

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

	/* problem size */
	n  = atoi(argv[1]);

	/* number of runs */
	num_of_runs = atoi(argv[2]);
	
	#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif		
	 
	if(rank == 0) wtimes = (Real*)malloc(sizeof(Real)*num_of_runs);

	for (i = 0; i < num_of_runs; ++i) {
		t1 = wall_time();
		umax = poisson_parallel(n, f ,u);	
		t2 = wall_time();
		if (rank == 0) {
			printf("Run: %d\t Time: %.16e\tUmax: %.16e\n", i, t2-t1, umax);
			wtimes[i] = t2- t1;
		}
	}
	

	/* find avarage walltime */
	if (rank == 0) {
 		printf("Average Time: %.16e\n", get_average(wtimes, num_of_runs, OUTLIERS_CUTOFF));		
	}
	
	#ifdef HAVE_MPI
	MPI_Finalize();
	#endif
	
	return 0;

}
