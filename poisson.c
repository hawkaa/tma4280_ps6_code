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
	
	MPI_Init(&argc, &argv);	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	 
	wtimes = (Real*)malloc(sizeof(Real)*num_of_runs);

	for (i = 0; i < num_of_runs; ++i) {
		umax = poisson_parallel(n, &wtimes[i], f ,u);	
		if (rank == 0) {
			printf("Run: %d\t Time: %.25e\tUmax: %.25e\n", i, wtimes[i], umax);
		}
	}
	

	/* find avarage walltime */
	if (rank == 0) {
 		printf("Average Time: %.25e\n", get_average(wtimes, num_of_runs,
						OUTLIERS_CUTOFF));		
	}
	
	return 0;

}
