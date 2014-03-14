/* global includes */
#include <stdlib.h>
#include <stdio.h>

/* local includes */
#include "ps6_common_library.h"


int
main(int argc, char** argv)
{
	int n;
	if (argc < 2  {
		printf("need a problem size\n");
		return 1;
	}
	n  = atoi(argv[1]);
	return poisson(n);

}
