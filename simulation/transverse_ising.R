#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(!is.na(args[1])){
    modus = args[[1]]
}else{
    modus = "plot"
}
source("ising_dos.R")

# physical parameters
J = 1
J2 = 0
h = seq(.5, 1.5, by=.1)

# auxiliary params
n.boot = 30
n_states = 49
steps = 1000000

# default lattice extent
t = 1i
L = 8
Nt = 10

# color coding for plots
cols = c("violet", "blue", "green", "orange", "red")

plot.dos = function(x, rho, var, var.name, obs, norm=FALSE, ...){
	scale = 1
	lapply(seq(rho), function(i){
			   if(norm) scale = exp(Re(rho[[i]]$norm))
			   plotwitherror(x=h, y=rho[[i]][[obs]][1,]*scale, dy=rho[[i]][[obs]][3,]*scale, mdy=rho[[i]][[obs]][2,]*scale, rep=i>1, col=cols[i], pch=i, xlab="h", ...)
	})
	legend("topright", col=cols, pch=1:5, legend=sapply(var, function(x) sprintf("%s=%g", var.name, x)))
}

plot.exact = function(h, var, FUN, ...){
	field = seq(min(h), max(h), length=1000)
	lapply(seq(var), function(i){
			   lines(field, sapply(field, function(x) FUN(h=x, ..., var[[i]])), col=cols[i])
	})
}

# vary L
L = 4 * 1:5
if(modus == "new"){
	res.L = lapply(L, function(x){
					 scan.w.err.op.midpoint(h=h, n.boot=n.boot, L=x, Nt=Nt, J=J, J2=J2, time=t, steps=steps, n_states=n_states)
	})
	save(res.L, file="res.L.RData")
}else if(modus == "plot" & file.exists("res.L.RData")){
	load("res.L.RData")
	plot.dos(h, res.L, L, "L", "real.sign", main="S_R of L", ylab="S_R", log="y", ylim=c(1e-4,1))
	plot.exact(h, L, free.real.sign, Nt=Nt, time=t)
	plot.dos(h, res.L, L, "L", "real.sign", norm=TRUE, main="free evolution of L", ylab="Tr(U)", log="y", ylim=c(1e-1,1e6))
	plot.exact(h, L, free.real.evolution, time=t)
	plot.dos(h, res.L, L, "L", "phase", main="S_I of L", ylab="S_I", log="y", ylim=c(1e-4,1))
	plot.dos(h, res.L, L, "L", "stat.pow", main="stat. pow. of L", ylab="stat. pow.", log="y", ylim=c(1e-4,1))
	plot.dos(h, res.L, L, "L", "stat.pow", norm=TRUE, main="Tr(U) of L", ylab="Tr(U)", log="y", ylim=c(1e-1,1e6))
}
L = 8

# vary t
t = .6i * 1:5
if(modus == "new"){
	res.t = lapply(t, function(x){
					 scan.w.err.op.midpoint(h=h, n.boot=n.boot, L=L, Nt=Nt, J=J, J2=J2, time=x, steps=steps, n_states=n_states)
	})
	save(res.t, file="res.t.RData")
}else if(modus == "plot" & file.exists("res.t.RData")){
	load("res.t.RData")
	plot.dos(h, res.t, Im(t), "t", "real.sign", main="S_R of t", ylab="S_R", log="y", ylim=c(1e-7,1))
	plot.exact(h, t, free.real.sign, Nt=Nt, L=L)
	plot.dos(h, res.t, Im(t), "t", "real.sign", norm=TRUE, main="free evolution of t", ylab="Tr(U)", log="y", ylim=c(1e-5,1e2))
	plot.exact(h, t, free.real.evolution, L=L)
	plot.dos(h, res.t, Im(t), "t", "phase", main="S_I of t", ylab="S_I", log="y", ylim=c(1e-2,1))
	plot.dos(h, res.t, Im(t), "t", "stat.pow", main="stat. pow. of t", ylab="stat. pow.", log="y", ylim=c(1e-3,1))
	plot.dos(h, res.t, Im(t), "t", "stat.pow", norm=TRUE, main="Tr(U) of t", ylab="Tr(U)", log="y", ylim=c(1e0,1e3))
}
t = 1i
L = 10

# vary Nt
Nt = 2 * 1:5
if(modus == "new"){
	res.Nt = lapply(Nt, function(x){
					 scan.w.err.op.midpoint(h=h, n.boot=n.boot, L=L, Nt=x, J=J, J2=J2, time=t, steps=steps, n_states=n_states)
	})
	save(res.Nt, file="res.Nt.RData")
}else if(modus == "plot" & file.exists("res.Nt.RData")){
	load("res.Nt.RData")
	plot.dos(h, res.Nt, Nt, "Nt", "real.sign", main="S_R of Nt", ylab="S_R", log="y", ylim=c(1e-4,1))
	plot.exact(h, Nt, free.real.sign, L=L, time=t)
	plot.dos(h, res.Nt, Nt, "Nt", "real.sign", norm=TRUE, main="free evolution of Nt", ylab="Tr(U)", log="y", ylim=c(1e-1,1e6))
	plot.exact(h, Nt, free.real.evolution, L=L, time=t)
	plot.dos(h, res.Nt, Nt, "Nt", "phase", main="S_I of Nt", ylab="S_I", log="y", ylim=c(1e-4,1))
	plot.dos(h, res.Nt, Nt, "Nt", "stat.pow", main="stat. pow. of Nt", ylab="stat. pow.", log="y", ylim=c(1e-4,1))
	plot.dos(h, res.Nt, Nt, "Nt", "stat.pow", norm=TRUE, main="Tr(U) of Nt", ylab="Tr(U)", log="y", ylim=c(1e-1,1.e6))
}
Nt = 10
