#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

typedef long long int ll;

double random_double(double minval, double maxval) {
	return minval + ((double)rand() / (double)RAND_MAX) * (maxval - minval);
}

int main (int argc, char *argv[])
{
	assert (argc == 2);
	ll N = atoll(argv[1]);
	srand(time(0));
	ll count = 0;
	ll times = N;
	while (times--) {
		double x = random_double(-1.,1);
		double y = random_double(-1.,1);
		if (x * x + y * y <= 1) count++;
	}
	printf("%.10lf\n", (double) count / (double) N * 4.);
}
