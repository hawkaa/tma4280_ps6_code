/* global includes */
#include <stdlib.h>
#include <stdio.h>

/* local includes */
#include "ps6_common_library.h"


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
	solution = poisson(n);
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			printf("%f\t", solution[i][j]);
		}
		printf("\n");
	}



}
