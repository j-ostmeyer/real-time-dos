#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "ising_aux.h"

double mod2pi(double x){
	// projects a phase into the interval [0:2pi]
	return x - M_2PI * floor(M_2PI_INV * x);
}

double phase2int(double x, unsigned n){
	// identifies [0:2pi] with the set {0...n}
	return n * (M_2PI_INV*x - floor(M_2PI_INV*x));
}

unsigned imag2int(double x, double max, unsigned n){
	return (unsigned)(n * (x/(2*max) + .5));
}

double add_log(double s1, double s2, double a, double b){
	// calculates exp(s) = a*exp(s1) + b*exp(s2)
	// returns s
	
	if(s2 < s1)
		return s1 + log(a + b*exp(s2-s1));
	else if(s1 < s2)
		return s2 + log(b + a*exp(s1-s2));
	else
		return s1 + log(a + b); // in case s1=s2=inf
}

unsigned strip_zero(double *x, unsigned n){
	unsigned i, count = 0;

	for(i = 0; i < n; i++)
		if(x[i] == 0){
			x[i] = -INFINITY;
			count++;
		}

	return count;
}

double normalise_log(double *x, unsigned n){
	unsigned i;
	double x_max = x[0], sum = 0;

	for(i = 1; i < n; i++) if(x[i] > x_max) x_max = x[i];

	for(i = 0; i < n; i++){
		x[i] -= x_max;
		sum += exp(x[i]);
	}
	sum = log(sum);

	for(i = 0; i < n; i++) x[i] -= sum;

	return x_max + sum;
}

void normalise_log2(double *x, unsigned n, double log_ratio){
	const double ratio = exp(log_ratio);
	const double a = ratio/(1+ratio);
	unsigned i;
	double norm;

	normalise_log(x, n);
	normalise_log(x+n, n);

	norm = log(a);
	for(i = 0; i < n; i++) x[i] += norm;

	n *= 2;
	norm = log(1-a);
	for(; i < n; i++) x[i] += norm;
}

double correlate(double *x, double *y, unsigned n){
	double sum = 0;
	unsigned i;
	for(i = 0; i < n; i++) sum += x[i]*y[i];
	return sum;
}

void print_vec(double *vec, unsigned n){
	unsigned i;
	for(i = 0; i < n; i++) printf("%g\n", vec[i]);
}

void print_cvec(double complex *vec, unsigned n){
	unsigned i;
	for(i = 0; i < n; i++) printf("%g\t%+g i\n", creal(vec[i]), cimag(vec[i]));
}

void print_mat(double *m, unsigned long n){
	unsigned i, k;

	for(i = 0; i < n; i++){
		for(k = 0; k < n; k++, m++) printf("%g,\t", *m);
		printf("\n");
	}
}

void print_spins(short *s, unsigned L, unsigned Nt){
	unsigned i, k;

	for(i = 0; i < L; i++){
		for(k = 0; k < Nt; k++) printf("%d\t", s[i*Nt+k]);
		printf("\n");
	}
}
