#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#define N 2000000
int nums[N];

double random_int(int minval, int maxval) {
	return minval + (rand() % (maxval - minval + 1));
}

int main (int argc, char *argv[])
{
	srand(time(0));
	long long int sum = 0;
	for (int i = 0; i < N; ++i) {
		nums[i] = random_int(-1000, 1000);
		sum += nums[i];
	}
	double mu = (double)sum / (double)N;
	double sigma = 0;
	for (int i = 0; i < N; ++i) {
		double tmp = nums[i] - mu;
		sigma += tmp * tmp;		
	}
	sigma = sqrt(sigma / N);
	printf("mu=%lf, sigma=%lf\n", mu, sigma);
}
