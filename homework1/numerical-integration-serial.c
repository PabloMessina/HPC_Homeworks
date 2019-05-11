#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#define N 2000000
typedef long long int ll;

int nums[N];

double PI, MU, SIGMA, _COEF1, _COEF2;

double random_int(int minval, int maxval) {
	return minval + (rand() % (maxval - minval + 1));
}

double random_double(double minval, double maxval) {
	return minval + ((double)rand() / (double)RAND_MAX) * (maxval - minval);
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
	for (int i = 2; i <= n - 2; i +=2) tmp += f(a + i*h);
	integral += 2*tmp;
	tmp = 0;
	for (int i = 1; i <= n - 1; i +=2) tmp += f(a + i*h);
	integral += 4*tmp;
	integral += f(b);

	return integral * h / 3.;
} 

int main (int argc, char *argv[])
{
		
	assert (argc == 3);
	srand(time(0));
	
	// ----- STEP 1: estimate PI constant -----
	ll count = 0;
	ll pi_iterations = atoll(argv[1]);
	for (ll i = 0; i < pi_iterations; ++i) {
		double x = random_double(-1.,1.);
		double y = random_double(-1.,1.);
		if (x * x + y * y <= 1) count++;
	}
	PI = (double)count / (double)pi_iterations * 4.;
	printf("PI %lf\n", PI);

	// ----- STEP 2: estimate MU and SIGMA -----
	MU = 0.;
	for (int i = 0; i < N; ++i) {
		nums[i] = random_int(-1000, 1000);
		MU += nums[i];
	}
	MU /= N;
	SIGMA = 0.;
	for (int i = 0; i < N; ++i) {
		double tmp = nums[i] - MU;
		SIGMA += tmp * tmp;
	}
	SIGMA = sqrt(SIGMA / N);
	printf("MU %lf\nSIGMA %lf\n", MU, SIGMA);

	// ----- STEP 3: compute numerical integration ----
	_COEF1 = 1. / sqrt(2. * PI * SIGMA * SIGMA);
	_COEF2 = - 1. / (2. * SIGMA * SIGMA);
	int steps = atoi(argv[2]);
	printf("integral %lf\n", _COEF1 * integrate(MU - SIGMA, MU + SIGMA, func, steps));
}
