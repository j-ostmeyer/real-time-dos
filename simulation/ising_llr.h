#ifndef ISING_H0
#define ISING_H0

int topology(short *s, unsigned *nnt, unsigned ns, unsigned nn);
int local_topology(short *s, unsigned *nnt, unsigned i, unsigned nn);
double complex global_h(short *s, unsigned *nnt, double complex *nnc, unsigned ns, unsigned nn);
double complex local_h(short *s, unsigned *nnt, double complex *nnc, unsigned i, unsigned nn);
void observable(short *s, double complex h, int topo, unsigned L, unsigned Nt, FILE *obs);

void sample_dos_imag(unsigned L, unsigned Nt, unsigned steps, unsigned thermalise, double complex *J, double complex J2, double complex h, int start, double *log_dos, int fix_topo, double log_topo, int avoid, unsigned n_states, double maxIm, double scale, double offset, int seed, FILE *obs);
#endif
