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
#define N_MAX 9

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
	int rank;
	#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#else
	rank = 0;
	#endif

	int i, j;
	Real umax;
	
	if (rank == 0) {
		printf("# Convergence test for Problem Set 6 poisson solver\n");
		printf("# Problem Size\tAbsolute error\n");
	}

	for (i = N_MIN; i <= N_MAX; ++i) {
		j = pow(2, i);
		umax = poisson(j, *f, *u);
		if (rank == 0) {
			printf("%i\t%.16e\n", j, umax);
		}
	}

	#ifdef HAVE_MPI
	MPI_Finalize();
	#endif
	
	exit(0);
}
