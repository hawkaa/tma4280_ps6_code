/* global includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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
	return 5 * pi * pi * sin(pi * x) * sin(2 * pi * y);
}

Real
u(Real x, Real y)
{
	Real pi;
	pi = 4.0 * atan(1.0);
	return sin(pi * x) * sin(2 * pi * y);
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
			x = (Real)(i) / (Real)(problem_size + 1);
			y = (Real)(j) / (Real)(problem_size + 1);
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
	int i, j;
	Real **b;
	Real umax;

	b = createReal2DArray(pow(2, N_MAX) - 1, pow(2, N_MAX) - 1);
	for (i = N_MIN; i <= N_MAX; ++i) {
		j = pow(2, i);
		b = poisson(j, *f);
		umax = get_umax(b, j, *u);
		printf("%i:\t%f\n", j, umax);
	}
	
	/*

	solution = poisson(n, *f);
	


	for (i = 0; i < n - 1; ++i) {
		for (j = 0; j < n - 1; ++j) {
			//printf("%f\t", solution[i][j]);
		}
		//printf("\n");
	}
	//printf("\n");
	Real x, y;
	Real umax, res;
	umax = 0.0;
	for (i = 1; i < n; ++i) {
		for (j = 1; j < n; ++j) {
			x = (Real)j / (Real)(n);
			y = (Real)i / (Real)(n);
			res = fabs(solution[i-1][j-1] - u(x, y));
			if ( res > umax) {
				umax = res;
			}
			//printf("%f\t", u(x, y));
		}
		//printf("\n");
	}
	printf("\nUMAX: %f\n", umax);
	*/

}
