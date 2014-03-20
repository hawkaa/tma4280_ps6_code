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
#define N_MAX 10

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

static Real
get_umax(Real **b, int problem_size, function2D reference)
{
	int i, j;
	Real umax, sum;
	Real x, y;
	umax = 0.0;
	for (i = 1; i < problem_size; ++i) {
		for (j = 1; j < problem_size; ++j) {
			x = (Real)(j) / (Real)(problem_size);
			y = (Real)(i) / (Real)(problem_size);
			sum = fabs((*reference)(x, y) - b[i-1][j-1]);
			if (sum > umax) {
				umax = sum;
			}
		}
	}
	return umax;

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

	printf("# Convergence test for Problem Set 6 poisson solver\n");
	printf("# Problem Size\tAbsolute error\n");

	for (i = N_MIN; i <= N_MAX; ++i) {
		j = pow(2, i);
		umax = poisson(j, *f);
		if (rank == 0) {
			printf("%i\t%.16e\n", j, umax);
		}
	}

	#ifdef HAVE_MPI
	MPI_Finalize();
	#endif
	
	exit(0);
}
