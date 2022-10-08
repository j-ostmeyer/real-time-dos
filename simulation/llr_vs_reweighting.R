#!/usr/bin/env Rscript
#options(warn=1)

source("ising_dos.R")

# auxiliary params
n.boot = 40
n_states = 10
steps = c(1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8)

# default lattice extent
L = LENGTH
times = ifelse(L <= 30, min(L,12), 3)
Nt = floor(4 * L / times)
t = 1i/times
avoid = AVOID

# physical parameters
J2 = .3
h = .6
dJ = CHAOS
J = read.table(dir(paste0("../production/L_", L, "_dJ_", dJ, "/"), pattern="*_field.txt", full.name=TRUE)[1])[,1]

t.scale = (dJ-1)*2/5 + 1
files = paste0("R.Data/DIR/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")

# vary t
k = OFFSET + 1:times

res.llr = lapply(steps, function(x){
					 scan.w.err.op.midpoint(fac=k, n.boot=n.boot, L=L, Nt=Nt, J=J, J2=J2, h=h, time=t/t.scale, steps=x, n_states=n_states, avoid=avoid)
	})
save(res.llr, file=files[1])

res.rew = lapply(steps, function(x){
					 scan.w.err.op.midpoint(fac=k, n.boot=n.boot, L=L, Nt=Nt, J=J, J2=J2, h=h, time=t/t.scale, steps=x, reweighting=TRUE, avoid=avoid)
	})
save(res.rew, file=files[2])
