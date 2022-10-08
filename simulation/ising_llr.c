#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#include "mt19937-64.h"
#include "ising_init.h"
#include "ising_aux.h"
#include "ising_lattice.h"
#include "ising_llr.h"

int topology(short *s, unsigned *nnt, unsigned ns, unsigned nn){
	// calculate number of "domain walls" in temporal direction
	// only works if number of lattice points divisible by 4
	unsigned i, pos;
	int top = 0;

	for(i = 0; i < ns; i++){
		pos = i*nn;
		top += s[i]*s[nnt[pos]];
	}

	return (ns-top)/4;
}

int local_topology(short *s, unsigned *nnt, unsigned i, unsigned nn){
	const unsigned pos = i*nn;
	const int env = s[nnt[pos]] + s[nnt[pos+1]];

	return -s[i]*env / 2;
}

double complex global_h(short *s, unsigned *nnt, double complex *nnc, unsigned ns, unsigned nn){
	// calculate the total Hamiltonian
	unsigned i, d, pos;
	double complex h = 0;

	for(i = 0; i < ns; i++){
		for(d = 0; d < nn; d += 2){ // don't double-count links
			pos = i*nn + d;
			h += s[i]*s[nnt[pos]] * nnc[pos];
		}
	}

	return h;
}

double complex local_h(short *s, unsigned *nnt, double complex *nnc, unsigned i, unsigned nn){
	// calculate the contribution to the Hamiltonian at position i
	unsigned d, pos;
	double complex h = 0;

	for(d = 0, pos = i*nn; d < nn; d++, pos++){
		h += s[nnt[pos]] * nnc[pos];
	}

	return h*s[i];
}

void observable(short *s, double complex h, int topo, unsigned L, unsigned Nt, FILE *obs){
	// dummy function set to whatever observable is required at compile time
	//
	// Hamiltonian
	fprintf(obs, "%.15g\t%.15g\t%d\n", creal(h), cimag(h), topo);
}

void sample_dos_imag(unsigned L, unsigned Nt, unsigned steps, unsigned thermalise, double complex *J, double complex J2, double complex field, int start, double *log_dos, int fix_topo, double log_topo, int avoid, unsigned n_states, double maxIm, double scale, double offset, int seed, FILE *obs){
	// given a lattice dimensionality and some couplings
	// calculate the densities of states of Im(H) with the probability density exp(-Re(H))
	const unsigned ns = L*Nt, nn = 6;
	unsigned long k;
	const unsigned long k_max = (unsigned long)steps * ns;
	double complex h, h_new;
	unsigned h_shift;
	unsigned h_ind, h_ind_new;
	int topo, topo_new;
	double weight = 0;

	unsigned *nnt = construct_lattice(L, Nt);
	double complex *nnc = construct_couplings(J, J2, field, L, Nt);
	init_genrand64(time(NULL)+seed);

	// initialisation
	short *s = init_lattice(ns, start);
	topo = topology(s, nnt, ns, nn);

	// thermalisation, i.e. updates disregarding Im(H)
	for(k = 0; k < thermalise*ns; k++){
		const unsigned j = (unsigned)(genrand64_real2()*ns);
		const double dH = -2*creal(local_h(s, nnt, nnc, j, nn));
		const int dTopo = -local_topology(s, nnt, j, nn);

		if((topo < avoid && dTopo > 0) || (topo+dTopo >= avoid && exp(dH) > genrand64_real2())){ // accept
			topo += dTopo;
			s[j] *= -1;
		}
	}

	h = global_h(s, nnt, nnc, ns, nn);
	h_shift = imag2int(cimag(h), maxIm, n_states);
	h_ind = h_shift + (1 & topo) * n_states;
	//printf("h = %g\n", h);

	// actual simulation
	for(k = 0; k < k_max; k++){
		if(k % ns == 0) observable(s, h, topo, L, Nt, obs);

		const unsigned j = (unsigned)(genrand64_real2()*ns);
		const double complex dH = -2*local_h(s, nnt, nnc, j, nn);
		const int dTopo = -local_topology(s, nnt, j, nn);

		h_new = h + dH;
		topo_new = topo + dTopo;
		h_shift = imag2int(cimag(h_new), maxIm, n_states);
		h_ind_new = h_shift + (1 & topo_new) * n_states;

		if(scale){
			weight = log_dos[h_ind] - log_dos[h_ind_new];
			weight -= ((1 & topo) - (1 & topo_new)) * log_topo;
		}

		if(topo_new >= avoid && genrand64_real2() < exp(creal(dH) + weight)){ // accept
			h = h_new;
			topo = topo_new;
			h_ind = h_ind_new;
			s[j] *= -1;
		}

		// update DoS or apply reweighting
		if(scale) log_dos[h_ind] += scale/(k+offset);
		else{
			const double complex phase = ((1&topo)?-1:1) * cexp(I*cimag(h));
			log_dos[0] += creal(phase);
			log_dos[1] += cimag(phase);
		}
	}

	// normalise DoS to 1
	if(scale){
		strip_zero(log_dos, 2*n_states);
		if(fix_topo) normalise_log2(log_dos, n_states, log_topo);
		else{
			for(k = 0; k < n_states; k++) log_dos[k] += log_topo;
			normalise_log(log_dos, 2*n_states);
		}
	}else{
		log_dos[0] /= k_max;
		log_dos[1] /= k_max;
	}

	free(s);
	free(nnt);
	free(nnc);
}
