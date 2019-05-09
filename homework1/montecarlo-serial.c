#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

double random_double(double minval, double maxval) {
	return minval + ((double)rand() / (double)RAND_MAX) * (maxval - minval);
}

int main (int argc, char *argv[])
{
	assert (argc == 2);
	int times = atoi(argv[1]);
	srand(time(0));
	int N = 0;
	int count = 0;
	while (times--) {
		double x = random_double(-1.,1);
		double y = random_double(-1.,1);
		if (x * x + y * y <= 1) count++;
		N++;
	}
	printf("%.4lf\n", (double) count / (double) N * 4.);
}
