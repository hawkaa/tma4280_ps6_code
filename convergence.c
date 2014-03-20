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
#define N_MIN 3
#define N_MAX 11	

/*
 * f function for the poisson problem
 * 5 * pi^2 * sin(pi * x) * sin(2 * pi * y)
 */
Real
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
Real
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
	/* get rank variable */
	int rank;
	#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#else
	rank = 0;
	#endif

	/* u_max variable */
	Real u_max;

	/* header data for gnuplot output */
	if (rank == 0) {
		printf("# Convergence test for Problem Set 6 poisson solver\n");
		printf("# Problem Size\tAbsolute error\n");
	}
	
	/* iterating over different problem sizes */
	for (i = N_MIN; i <= N_MAX; ++i) {
		j = pow(2, i);
		umax = poisson_parallel(j, *f, *u);
		if (rank == 0) {
			/* only rank 0 have valid result, and will print to file */
			printf("%i\t%.16e\n", j, umax);
		}
	}

	/* finalize mpi */
	#ifdef HAVE_MPI
	MPI_Finalize();
	#endif
	
	exit(0);
}
