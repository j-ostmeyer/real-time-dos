dyn.load("ising_dos.so")
source("ising_analytic.R")
source("ising_iterate_dos.R")
source("ising_obs.R")
source("ising_spline.R")

## Simulate the classical L x Nt Ising model using LLR or reweighting
sample_dos = function(L, Nt, steps, thermalise=ceiling(2*log(L*Nt)), J, J2=0, h=J, start="cold", n_states=L*Nt%/%2+1, log_dos=rep(0,2*n_states), fix.topo=FALSE, log.topo, avoid=0, offset=3, estimate.occupied=1, scale=.5*L*Nt*n_states*estimate.occupied/log(steps/offset), seed=0, observable="/dev/null", ...){
	stopifnot(length(J) == 1 || length(J) == L)
	stopifnot(length(h) == 1)
	stopifnot(length(log_dos) == 2*n_states)
	if(length(J) == 1){
		J = rep(J, L)
	}

	if(start == "cold") st = 1
	else if(start == "hot") st = 0
	else stop("Start type not supported!")

	if(missing(log.topo)) log.topo = 0
	else if(missing(fix.topo)) fix.topo = TRUE

	if(!missing(seed) & !missing(observable)) observable = paste0("s_", seed, "_", observable)

	max.im = Nt * sum(abs(J) + abs(J2))
	max.im = max.im * (1 + 1/(n_states-1)) # make physical maximum center of last bin
	im.S = (1:n_states - .5)*2*max.im/n_states - max.im

	dos = .Call("sample_dos", L, Nt, as.integer(steps), thermalise, as.complex(J), as.complex(J2), as.complex(h), st, as.vector(log_dos), as.integer(fix.topo), log.topo, avoid, n_states, max.im, scale, offset*n_states, seed, observable)

	return(list(ImS = im.S, dos = matrix(exp(dos), ncol=2), steps = steps, offset = offset+steps*L*Nt/n_states, scale = scale, observable = observable))
}

## Simulate the quantum length-L transverse Ising model using LLR or reweighting
transverse_ising_dos = function(L, Nt, J, J2=0, h, time, log.topo.ratio, avoid=0, ana.offset, reweighting=FALSE, observable="/dev/null", ...){
	delta = time/Nt
	h_cl = -1/2*log(tanh(delta*h))
	if(Re(delta*h) == 0){ # pure real time evolution, take care of sign by topology
		h_cl = Re(h_cl)
	}

	sign.problem = free.real.sign(L, Nt, h, time)
	norm = exp(Re(norm.dos(L, Nt, h, time, quantum=FALSE)))
	if(missing(log.topo.ratio)){
		if(avoid == 0) log.topo.ratio = log((1+sign.problem) / (1-sign.problem))
		if(avoid == 1) log.topo.ratio = log((1+sign.problem - 2*2^L/norm) / (1-sign.problem))
		if(avoid == 2) log.topo.ratio = log((1+sign.problem - 2*2^L/norm) / (1-sign.problem - 2*2^L*tan(abs(delta*h))^2*Nt*(Nt-1)/2*L/norm))
	}

	if(reweighting){
		res = sample_dos(L, Nt, J=delta*J, J2=delta*J2, h=h_cl, n_states=1, scale=0, avoid=avoid, ...)
		res$dos = log(res$dos[1,])
	}else{
		res = sample_dos(L, Nt, J=delta*J, J2=delta*J2, h=h_cl, log.topo=log.topo.ratio, avoid=avoid, ...)
	}

	res$L = L
	res$Nt = Nt
	res$J = J
	res$J2 = J2
	res$h = h
	res$h_cl = h_cl
	res$time = time

	res$avoid = avoid
	if(avoid == 0) res$renorm = 1
	if(avoid == 1) res$renorm = 1 - 2^L/norm
	if(avoid == 2) res$renorm = 1 - 2^L*(1 + tan(abs(delta*h))^2*Nt*(Nt-1)/2*L)/norm
	res$dos = res$dos * res$renorm

	res$sign.problem = sign.problem
	res$topo.ratio = exp(log.topo.ratio)
	res$norm = exp(Re(norm.dos(L, Nt, h, time)))
	res$c.norm = norm

	if(missing(ana.offset)) res$ana.offset = leading.order.sff(res, J0=1)
	else res$ana.offset = ana.offset/res$c.norm

	if(reweighting){
		if(missing(observable)){
			st.p = res$ana.offset + res$dos[1] + res$dos[2]*1i
			res$meas = c(NA, NA, st.p, st.p*res$norm, NA, NA)
		}else{
			energies = read.table(observable)
			res$meas = order.params.reweighting(energies, res)
			file.remove(observable)
		}
	}else{
		res$meas = order.params.midpoint(res)
	}

	return(res)
}
