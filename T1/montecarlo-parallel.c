#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

double random_double(double minval, double maxval, unsigned int *seed) {
	return minval + ((double)rand_r(seed) / (double)RAND_MAX) * (maxval - minval);
}

int main (int argc, char *argv[])
{
	// make sure we receive 2 arguments
	assert (argc == 3);
	long long int N = atoll(argv[1]); // N = number of iterations
	int n_threads = atoi(argv[2]); // number of threads
	// printf("N = %lld, n_threads = %d\n", N, n_threads); // debugging

	// init different seeds for each thread
	srand(time(0));
	int seeds[n_threads];
	for (int i = 0; i < n_threads; ++i) seeds[i] = rand();

	int count = 0; // global counter
	omp_set_dynamic(0);
	omp_set_num_threads(n_threads);
	#pragma omp parallel default(none) shared(count, seeds, N, n_threads)
	{
		int thread_i = omp_get_thread_num(); // get thread number
		assert (thread_i < n_threads); 
		unsigned int s = seeds[thread_i]; // obtain local seed
		
		#pragma omp for reduction(+:count)
		for (int i = 0; i < N; ++i) {
			double x = random_double(-1.,1, &s);
			double y = random_double(-1.,1, &s);
			if (x * x + y * y <= 1) count++;
		}
	}
	printf("%.10lf\n", (double)count / (double) N * 4.);
}
