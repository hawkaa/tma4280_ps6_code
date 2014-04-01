/* global includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/* local includes */
#include "ps6_common_library.h"

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
	/* u_max variable */
	Real u_max;

	/* loop variables */
	int j, i;

	/* read arguments */
	if (argc < 3) {
		printf("Need to input n_min and n_max for convergence test\n");
		exit(1);
	}
	int n_min = atoi(argv[1]);
	int n_max = atoi(argv[2]);

	/* get rank variable */
	int rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);




	/* header data for gnuplot output */
	if (rank == 0) {
		printf("# Convergence test for Problem Set 6 poisson solver\n");
		printf("# Problem Size\tAbsolute error\n");
	}
	
	/* iterating over different problem sizes */
	for (i = n_min; i <= n_max; ++i) {
		j = pow(2, i);
		u_max = poisson_parallel(j, NULL, *f, *u);
		if (rank == 0) {
			/* only rank 0 have valid result, and will print to file */
			printf("%i\t%.25e\n", j, u_max);
		}
	}

	/* finalize mpi */
	MPI_Finalize();
	
	exit(0);
}
