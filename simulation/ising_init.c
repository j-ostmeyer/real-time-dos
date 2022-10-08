#include <stdio.h>
#include <stdlib.h>

#include "mt19937-64.h"
#include "ising_init.h"

short *init_lattice(unsigned ns, int start){
	unsigned i;
	short *s = malloc(ns*sizeof(short));

	if(start == 1){ // cold start
		for(i = 0; i < ns; i++) s[i] = 1;
	}else{ // hot start
		for(i = 0; i < ns; i++) s[i] = (genrand64_real2() > 0.5)? 1:-1;
	}

	return s;
}
