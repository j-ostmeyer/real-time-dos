#include <R.h>
#include <Rinternals.h>
#include <complex.h>

#include "ising_llr.h"
#include "ising_noisy_trace.h"
#include "ising_analytic.h"

static R_INLINE double complex toC99(const Rcomplex *x){
#if __GNUC__
	double complex ans = (double complex) 0; /* -Wall */
	__real__ ans = x->r;
	__imag__ ans = x->i;
	return ans;
#else
	return x->r + x->i * I;
#endif
}

static R_INLINE void SET_C99_COMPLEX(Rcomplex *x, R_xlen_t i, double complex value){
	Rcomplex *ans = x+i;
	ans->r = creal(value);
	ans->i = cimag(value);
}

SEXP sample_dos(SEXP L, SEXP Nt, SEXP steps, SEXP thermalise, SEXP J, SEXP J2, SEXP h, SEXP start, SEXP ext_log_dos, SEXP fix_topo, SEXP log_topo, SEXP avoid, SEXP n_states, SEXP maxIm, SEXP scale, SEXP offset, SEXP seed, SEXP observable){
	const unsigned n = asInteger(n_states);
	SEXP out = PROTECT(allocVector(REALSXP, 2*n)); // need space for rho

	FILE *obs = fopen(CHAR(asChar(observable)), "w");

	double complex *field = malloc(asInteger(L)*sizeof(double complex));
	for(int i = 0; i < asInteger(L); i++) field[i] = toC99(COMPLEX(J) + i);

	double *log_dos = calloc(2*n, sizeof(double));
	if(LENGTH(ext_log_dos)){
		memcpy(log_dos, REAL(ext_log_dos), 2*n*sizeof(double));
	}

	sample_dos_imag(asInteger(L), asInteger(Nt), asInteger(steps), asInteger(thermalise), field, toC99(COMPLEX(J2)), toC99(COMPLEX(h)), asInteger(start), log_dos, asInteger(fix_topo), asReal(log_topo), asInteger(avoid), n, asReal(maxIm), asReal(scale), asReal(offset), asInteger(seed), obs);

	memcpy(REAL(out), log_dos, 2*n*sizeof(double));

	free(field);
	free(log_dos);
	fclose(obs);
	UNPROTECT(1);

	return out;
}

SEXP noisy_trace(SEXP L, SEXP J, SEXP h, SEXP t, SEXP t_step, SEXP n_sources, SEXP scheme, SEXP time_series, SEXP full_data, SEXP operator_id, SEXP first_all_x, SEXP seed){
	const unsigned long n_steps = (unsigned long) ceil(cabs(toC99(COMPLEX(t))) / asReal(t_step));
	const unsigned long n = (asInteger(time_series)? n_steps:1L) * (asInteger(full_data)? asInteger(n_sources):2);
	SEXP out = PROTECT(allocVector(CPLXSXP, n)); // need space for rho

	double *coupling = malloc(3*sizeof(double));
	memcpy(coupling, REAL(J), 3*sizeof(double));
	double *field = malloc(asInteger(L)*sizeof(double));
	memcpy(field, REAL(h), asInteger(L)*sizeof(double));

	double complex *trace;
	trace = trace_estimator(coupling, field, toC99(COMPLEX(t)), asInteger(L), n_steps, asInteger(n_sources), asInteger(scheme), asInteger(time_series), asInteger(full_data), asInteger(operator_id), asInteger(first_all_x), asInteger(seed));

	for(unsigned long i = 0; i < n; i++) SET_C99_COMPLEX(COMPLEX(out), i, trace[i]);

	free(coupling);
	free(field);
	free(trace);
	UNPROTECT(1);

	return out;
}

SEXP classical_approx(SEXP L, SEXP Nt, SEXP J, SEXP J2, SEXP h, SEXP J0, SEXP t, SEXP avoid){
	SEXP out = PROTECT(allocVector(CPLXSXP, 1)); // need space for rho
	double *field = malloc(asInteger(L)*sizeof(double));
	memcpy(field, REAL(J), asInteger(L)*sizeof(double));

	double complex res;
	res = leading_order(asInteger(L), asInteger(Nt), field, asReal(J2), asReal(h), asReal(J0), toC99(COMPLEX(t)), asInteger(avoid));

	SET_C99_COMPLEX(COMPLEX(out), 0, res);

	free(field);
	UNPROTECT(1);

	return out;
}
