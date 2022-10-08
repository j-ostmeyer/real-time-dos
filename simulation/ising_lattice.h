#ifndef ISING_LATTICE
#define ISING_LATTICE

unsigned *construct_lattice(unsigned L, unsigned Nt);
double complex *construct_couplings(double complex *J, double complex J2, double complex h, unsigned L, unsigned Nt);
#endif
