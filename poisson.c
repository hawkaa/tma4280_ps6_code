/* global includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

int
main(int argc, char** argv)
{
	int n, i, j;
	Real **solution;
	if (argc < 2)  {
		printf("need a problem size\n");
		return 1;
	}
	n  = atoi(argv[1]);
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


}
