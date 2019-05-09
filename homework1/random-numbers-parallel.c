#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#define N 2000000
int nums[N];

double random_int(int minval, int maxval, int *seed) {
	return minval + (rand_r(seed) % (maxval - minval + 1));
}

int main (int argc, char *argv[])
{
	assert (argc == 2);
	int n_threads = atoi(argv[1]);

	int seeds[n_threads];
	srand(time(0));
	for (int i = 0; i < n_threads; ++i) seeds[i] = rand();
	long long int sum = 0;

	omp_set_dynamic(0);
	omp_set_num_threads(n_threads);	
	#pragma omp parallel default(none) shared(nums, sum, seeds, n_threads)
	{
		int thread_i = omp_get_thread_num();
		assert (thread_i < n_threads);
		unsigned int s = seeds[thread_i];

		#pragma omp for reduction(+:sum) schedule(static)
		for (int i = 0; i < N; ++i) {
			nums[i] = random_int(-1000, 1000, &s);
			sum += nums[i];
		}
	}
	
	double mu = (double)sum / (double)N;
	double sigma = 0;
	
	#pragma omp parallel default(none) shared(nums, mu, sigma)
	{
		#pragma omp for reduction(+:sigma) schedule(static)
		for (int i = 0; i < N; ++i) {
			double tmp = nums[i] - mu;
			sigma += tmp * tmp;
		}
	}

	sigma = sqrt(sigma / N);
	printf("mu=%lf, sigma=%lf\n", mu, sigma);
}
