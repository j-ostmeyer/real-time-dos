
get.ImS = function(J, J2, time, n_states){
	max.im = abs(time) * sum(abs(J) + abs(J2))
	max.im = max.im * (1 + 1/(n_states-1)) # make physical maximum center of last bin
	im.S = (1:n_states - .5)*2*max.im/n_states - max.im
	return(im.S)
}

sd.midpoint = function(rho){
	s.sq = rho$ImS^2
	return(apply(rho$dos, 2, function(x) sum(x*s.sq)/sum(x)))
}

stat.power.midpoint = function(rho, J0=1, range){
	s.exp = exp(1i*rho$ImS[range]*J0)
	phase = apply(rho$dos[range,], 2, function(x) sum(s.exp*x))
	return(c(phase, phase[1]-phase[2]))
}

order.params.midpoint = function(rho, J0=1, cut.off){
	n = dim(rho$dos)[1]
	if(missing(cut.off)) ii = 1:n
	else if(cut.off < 1) ii = which(abs((rho$dos[,1]-rho$dos[,2])/(rho$dos[,1]+rho$dos[,2])) > cut.off)
	else ii = unique(c(1:cut.off, (n+1-cut.off):n))

	real.sign = sum(rho$dos[ii,1] - rho$dos[ii,2])
	st.p = stat.power.midpoint(rho, J0, ii)
	phase = sum(st.p[1:2])

	if(is.null(rho$ana.offset) | J0 != 1)
		ana.offset = leading.order.sff(rho, J0)
	else ana.offset = rho$ana.offset
	st.p[3] = st.p[3] + ana.offset

	if(is.null(rho$hermit.fit))
		hermit.res = NA
	else
		hermit.res = hermit.ft(rho$hermit.fit$t0)

	return(c(real.sign, phase, st.p[3], st.p[3]*rho$norm, hermit.res, (hermit.res+ana.offset)*rho$norm))
}

order.params.reweighting = function(energies, rho, J0=1){
	real.sign = mean(ifelse(energies[,3]%%2, -1, 1))
	phase = mean(exp(1i*J0*energies[,2]))
	st.p = mean(exp(1i*J0*energies[,2])*ifelse(energies[,3]%%2, -1, 1))*rho$renorm
	st.p = st.p + leading.order.sff(rho, J0)

	return(c(real.sign, phase, st.p, st.p*rho$norm, NA, NA))
}

scan.op.midpoint = function(h, J0=1, ...){
	res = sapply(h, function(x){
					 rho = transverse_ising_dos(..., h=x)
					 out = order.params.midpoint(rho, J0)
					 return(out)
				})
	return(res)
}

scan.w.err.op.midpoint = function(time, Nt, n_states, fac, n.boot=10, ...){
	reweighting = missing(n_states)
	measure = function(x, seed, ...){
		if(reweighting)
			rho = transverse_ising_dos(..., time=time*x, Nt=Nt*x, seed=seed)
		else
			rho = transverse_ising_dos(..., time=time*x, Nt=Nt*x, n_states=n_states*x, seed=seed)
		return(c(Re(rho$meas[4]), Im(rho$meas[4]), as.vector(rho$dos)))
	}
	res = lapply(fac, function(x){
					 ana.offset = classical_approx(..., Nt=Nt*x, t=time*x)
					 bootstrap.asym(1000*(1:n.boot), measure, x=x, ana.offset=ana.offset, ...)
				})

	re = sapply(res, function(x) x$res[,1])
	im = sapply(res, function(x) x$res[,2])

	sff = sapply(res, function(r){
					 x = r$res
					 c(sqrt(x[1,1]^2 + x[1,2]^2), sqrt(((x[1,1]*x[2:3,1])^2 + (x[1,2]*x[2:3,2])^2)/(x[1,1]^2 + x[1,2]^2)))
				})

	dos = lapply(res, function(x) x$res[,-(1:2)])

	meas.samples = lapply(res, function(x) x$samples[1:2,])
	dos.samples = lapply(res, function(x) x$samples[-(1:2),])

	scan = list(re=re, im=im, sff=sff, dos=dos, meas.samples=meas.samples, dos.samples=dos.samples)
	return(scan)
}

extend.scan = function(scan, fac, J, J2, time, n_states){
	spline = lapply(seq(fac), function(i){
						spline.sff.w.err(get.ImS(J, J2, time*fac[i], n_states*fac[i]), scan$dos[[i]], scan$dos.samples[[i]])
					})
	scan$spline = spline

	return(scan)
}


scan.w.err.dos = function(time, Nt, n_states, fac, n.boot=10, ...){
	measure = function(x, seed, ...){
		rho = transverse_ising_dos(..., time=time*x, Nt=Nt*x, n_states=n_states*x, seed=seed)
		return(as.vector(rho$dos))
	}
	res = lapply(fac, function(x){
					 bootstrap.asym(1000*(1:n.boot), measure, x=x, ...)
				})

	return(res)
}

# An mc-version of the sapply function.
mcsapply <- function (X, FUN, ..., simplify=TRUE, USE.NAMES=TRUE,
					  mc.cores=min(length(X), parallel::detectCores()), mc.preschedule=FALSE,
					  affinity.list=seq(X)%%mc.cores+1) {
	FUN <- match.fun(FUN)
	answer <- parallel::mclapply(X = X, FUN = FUN, ...,
								 mc.cores=mc.cores, mc.preschedule=mc.preschedule, affinity.list=affinity.list)
	if (USE.NAMES && is.character(X) && is.null(names(answer))) 
		names(answer) <- X
	if (!isFALSE(simplify) && length(answer)) 
		simplify2array(answer, higher = (simplify == "array"))
	else answer
}

bootstrap.asym = function(seed, FUN, ...){
	# FUN takes a seed and some other arguments and returns a scalar or vector
	# for all results of FUN the median and sigma-quantiles are returned
	#             median ...
	# 3xn-matrix: q(16%) ...
	#             q(84%) ...
	samples = mcsapply(seed, function(s, ...){
						 FUN(..., seed=s)
				}, ...)
	if(is.null(dim(samples))){
		samples = t(samples)
	}

	res = apply(samples, 1, quantile, probs=c(.5, .16, .84), names=FALSE, na.rm=TRUE)
	res[2,] = (res[1,] - res[2,])/sqrt(length(seed)) # negative error estimation
	res[3,] = (res[3,] - res[1,])/sqrt(length(seed)) # positive error estimation

	return(list(res=res, samples=samples))
}

log.smoothen = function(x, b=length(x), n=2*log(b)/log(2)){
	y = sapply(seq(n), function(i) mean(x[floor(b^((i-1)/n)):floor(b^(i/n))]))
	return(y)
}
