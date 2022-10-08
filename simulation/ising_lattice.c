#include <stdlib.h>
#include <complex.h>
#include "ising_lattice.h"

unsigned *construct_lattice(unsigned L, unsigned Nt){
	// constructs nearest neighbour table for rectangle with next-to-nearest neighbour coupling
	// 1. lattice has time extend Nt
	// 2. spin chain has length L
	// 3. two additional nn neighbours in spatial direction
	const unsigned dim = 2;
	const unsigned size = L*Nt;
	const unsigned neighbours = 2*(dim+1);
	unsigned d, i;
	unsigned l[3], N[3];
	unsigned *nnt;
	int id;

	l[0] = Nt;
	l[1] = L;

	N[0] = 1;
	for(d = 1; d < dim; d++){
		N[d] = N[d-1]*l[d-1];
	}

	l[d] = l[d-1];
	N[d] = N[d-1];

	nnt = malloc(neighbours * size * sizeof(unsigned));

	for(i = 0; i < size; i++){
		for(d = 0; d < dim; d++){
			id = (i/N[d]) % l[d];
			nnt[i*neighbours + 2*d] = i + ((id+1)%l[d] - id)*N[d];
			nnt[i*neighbours + 2*d + 1] = i + ((id+l[d]-1)%l[d] - id)*N[d];
		}
		nnt[i*neighbours + 2*d] = i + ((id+2)%l[d] - id)*N[d];
		nnt[i*neighbours + 2*d + 1] = i + ((id+l[d]-2)%l[d] - id)*N[d];
	}

	return nnt;
}

double complex *construct_couplings(double complex *J, double complex J2, double complex h, unsigned L, unsigned Nt){
	// constructs coupling table
	// 1. h is constant coupling in temporal direction
	// 2. J is variable coupling in spatial direction
	// 3. J2 is constant nn neighbour coupling in spatial direction
	const unsigned dim = 2;
	const unsigned size = L*Nt;
	const unsigned neighbours = 2*(dim+1);
	unsigned i;
	double complex *nnc;
	int id;

	nnc = malloc(neighbours * size * sizeof(double complex));

	for(i = 0; i < size; i++){
		nnc[i*neighbours] = h;
		nnc[i*neighbours + 1] = h;

		id = i/Nt;
		nnc[i*neighbours + 2] = J[id];
		nnc[i*neighbours + 3] = J[(id+L-1)%L];

		nnc[i*neighbours + 4] = J2;
		nnc[i*neighbours + 5] = J2;
	}

	return nnc;
}
