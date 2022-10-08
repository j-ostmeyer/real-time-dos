#ifndef ISING_AUX
#define ISING_AUX

#define M_2PI 6.283185307179586
#define M_2PI_INV 0.1591549430918953

double mod2pi(double x);
double phase2int(double x, unsigned n);
unsigned imag2int(double x, double max, unsigned n);

double add_log(double s1, double s2, double a, double b);
unsigned strip_zero(double *x, unsigned n);
double normalise_log(double *x, unsigned n);
void normalise_log2(double *x, unsigned n, double log_ratio);

double correlate(double *x, double *y, unsigned n);

void print_vec(double *vec, unsigned n);
void print_cvec(double complex *vec, unsigned n);
void print_mat(double *m, unsigned long n);
void print_spins(short *s, unsigned L, unsigned Nt);

#endif
