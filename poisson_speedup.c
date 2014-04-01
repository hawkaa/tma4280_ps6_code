/* global includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/* local includes */
#include "ps6_common_library.h"

/* constants */
#define OUTLIERS_CUTOFF 2

/*
 * f function for the poisson problem
 * 5 * pi^2 * sin(pi * x) * sin(2 * pi * y)
 */
static Real
f(Real x, Real y)
{
	Real pi;
	pi = 4.0 * atan(1.0);
	return 5.0 * pi * pi * sin(pi * x) * sin(2.0 * pi * y);
}

/*
 * u (reference) function for the poisson problem
 * 5 * pi^2 * sin(pi * x) * sin(2 * pi * y)
 */
static Real
u(Real x, Real y)
{
	Real pi;
	pi = 4.0 * atan(1.0);
	return sin(pi * x) * sin(2.0 * pi * y);
}

/*
 * Main function 
 */
int
main(int argc, char** argv)
{
	/* loop variables */
	int j, i, p;

	/* read arguments */
	if (argc < 4) {
		printf("Need to input n_min and n_max for speedup test\n");
		exit(1);
	}

	int num_runs = atoi(argv[1]);
	int n_min = atoi(argv[2]);
	int n_max = atoi(argv[3]);

	/* get rank variable */
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	Real *wtimes = (Real*)malloc(sizeof(Real) * num_runs);

	/* header data for gnuplot output */
	if (rank == 0) {
		printf("# Speedup test for Problem Set 6 poisson solver\n");
		printf("# Problem Size\tTime\n");
	}
	
	/* iterating over different problem sizes */
	for (i = n_min; i <= n_max; ++i) {
		j = pow(2, i);
		for (p = 0; p < num_runs; ++p) {
			poisson_parallel(j, &wtimes[p], *f, NULL);
		}
		if (rank == 0) {
			/* only rank 0 have valid result, and will print to 
			file */
			printf("%i\t%.25e\n", j, get_average(wtimes, num_runs,
			OUTLIERS_CUTOFF));
		}
	}

	/* finalize mpi */
	MPI_Finalize();
	
	exit(0);
}
