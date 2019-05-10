#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#define N 2000000
int nums[N];

double PI, MU, SIGMA, _COEF1, _COEF2;

double random_int(int minval, int maxval, int* seed) {
	return minval + (rand_r(seed) % (maxval - minval + 1));
}

double random_double(double minval, double maxval, int* seed) {
	return minval + ((double)rand_r(seed) / (double)RAND_MAX) * (maxval - minval);
}

double func(double x) {
	double tmp = x - MU;
	return exp(_COEF2 * tmp * tmp);
}

double integrate(double a, double b, double (*f)(double), int n) {
	double h = (b - a)/n;
	assert (h > 0);

	double integral = 0, tmp;
	integral += f(a);
	tmp = 0;
	#pragma omp parallel default(none) shared(n, tmp, f, a, h)
	{
		#pragma omp for reduction(+:tmp)
		for (int i = 2; i <= n - 2; i +=2) tmp += f(a + i*h);
	}
	integral += 2*tmp;
	tmp = 0;
	#pragma omp parallel default(none) shared(n, tmp, f, a, h)
	{
		#pragma omp for reduction(+:tmp)
		for (int i = 1; i <= n - 1; i +=2) tmp += f(a + i*h);
	}
	integral += 4*tmp;
	integral += f(b);

	return integral * h / 3.;
} 

int main (int argc, char *argv[])
{
		
	assert (argc == 4);
	srand(time(0));
	int n_threads = atoi(argv[3]);
	omp_set_dynamic(0);
	omp_set_num_threads(n_threads);
	
	// ----- STEP 0: define different seeds per thread ---
	srand(time(0));
	int seeds[n_threads];
	for (int i = 0; i < n_threads; ++i) seeds[i] = rand();
	
	// ----- STEP 1: estimate PI constant -----
	int count = 0;	
	int pi_iterations = atoi(argv[1]);
	#pragma omp parallel default(none) shared(pi_iterations, count, seeds)
	{
		unsigned int s = seeds[omp_get_thread_num()];
		#pragma omp for reduction(+:count)
		for (int i = 0; i < pi_iterations; ++i) {
			double x = random_double(-1.,1, &s);
			double y = random_double(-1.,1, &s);
			if (x * x + y * y <= 1) count++;
		}
	}
	PI = (double)count / (double)pi_iterations * 4.;
	printf("PI = %lf\n", PI);

	// ----- STEP 2: estimate MU and SIGMA -----
	MU = 0.;
	#pragma omp parallel default(none) shared(MU, nums, seeds)
	{
		unsigned int s = seeds[omp_get_thread_num()];
		#pragma omp for reduction(+:MU) schedule(static)
		for (int i = 0; i < N; ++i) {
			nums[i] = random_int(-1000, 1000, &s);
			MU += nums[i];
		}
	}
	MU /= N;
	SIGMA = 0.;
	#pragma omp parallel default(none) shared(SIGMA, nums, MU)
	{
		#pragma omp for reduction(+:SIGMA) schedule(static)
		for (int i = 0; i < N; ++i) {
			double tmp = nums[i] - MU;
			SIGMA += tmp * tmp;
		}
	}
	SIGMA = sqrt(SIGMA / N);
	printf("MU=%lf, SIGMA=%lf\n", MU, SIGMA);

	// ----- STEP 3: computer numerica integration ----
	_COEF1 = 1. / sqrt(2. * PI * SIGMA * SIGMA);
	_COEF2 = - 1. / (2. * SIGMA * SIGMA);
	int steps = atoi(argv[2]);
	printf("integral = %lf\n", _COEF1 * integrate(MU - SIGMA, MU + SIGMA, func, steps));
}
