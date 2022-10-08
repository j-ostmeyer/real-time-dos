#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

double sinc(double x){
	return x? sin(x)/x : 1;
}

double sinsin(double x, double a, double b){
	return x? sin(a*x)/sin(b*x) : a/b;
}

double complex leading_order(unsigned L, unsigned Nt, double *J, double J2, double h, double J0, double complex time, unsigned avoid){
	const unsigned long N = 1L << L;
	if(time == 0) return N;

	const double delta = cabs(time)/Nt, delta2 = 2*delta;
	const double facT = .5 * (Nt? pow(tan(delta*h), 2) * Nt: pow(cabs(time)*h, 2));
	const double facS = delta2 * (Nt-1);
	const double time2 = 2*cabs(time);
	const double complex t = J0*time, t2 = -2*t;

	double resR = 0, resI = 0;

#pragma omp parallel for reduction(+:resR,resI)
	for(unsigned long k = 0; k < N; k++){
		double energy = 0;
		double complex phase = 0;

		for(unsigned i = 0; i < L; i++){
			unsigned sign;
			double coupling = 0;

			sign = ((k >> i) ^ (k >> ((i+1)%L))) & 1;
			coupling += sign? -J[i]:J[i];
			sign = ((k >> i) ^ (k >> ((i+L-1)%L))) & 1;
			coupling += sign? -J[(i+L-1)%L]:J[(i+L-1)%L];

			sign = ((k >> i) ^ (k >> ((i+2)%L))) & 1;
			coupling += sign? -J2:J2;
			sign = ((k >> i) ^ (k >> ((i+L-2)%L))) & 1;
			coupling += sign? -J2:J2;

			coupling *= .5;

			energy += coupling;
			if(avoid == 2)
				phase += cexp(t2*coupling) * (Nt? sinsin(coupling, facS, delta2) : sinc(time2*coupling));
		}

		const double complex res = cexp(t*energy) * (1 - facT*phase);
		resR += creal(res);
		resI += cimag(res);
	}

	return resR + I*resI;
}
