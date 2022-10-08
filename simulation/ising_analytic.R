
brute.force = function(L, J, J2=0, split=FALSE){
	div = 2^(0:(L-1))
	e = sapply(0:(2^L-1), function(i){
				   s = ifelse((i %/% div) %% 2, 1, -1)
				   s1 = c(s[2:L], s[1])
				   s2 = c(s[3:L], s[1:2])
				   if(split){
					   if(length(J) == L) Jm1 = c(J[L], J[1:(L-1)])
					   else Jm1 = J
					   if(length(J2) == L) J2m2 = c(J2[(L-1):L], J2[1:(L-2)])
					   else J2m2 = J2
					   sm1 = c(s[L], s[1:(L-1)])
					   sm2 = c(s[(L-1):L], s[1:(L-2)])
					   return((J*s1 + Jm1*sm1 + J2*s2 + J2m2*sm2) * s/2)
				   }else return(sum(J*s*s1 + J2*s*s2))
				})
	return(e)
}

classical_approx = function(L, Nt, J, J2, h, t, J0=1, avoid=2, ...){
	if(avoid == 0) return(0)
	return(.Call("classical_approx", L, Nt, J, J2, h, J0, as.complex(t), avoid))
}

leading.order.sff = function(rho, J0){
	lo = classical_approx(rho$L, rho$Nt, rho$J, rho$J2, rho$h, rho$time, J0, rho$avoid)
	return(lo/rho$c.norm)
}

log.z.classical.1d = function(h, Nt, boundaries="periodic"){
	c = log(cosh(h))
	t = tanh(h)
	switch(boundaries,
		   "periodic" = return(Nt*log(2) + Nt*c + log(1+t^Nt)),
		   "open" = return(Nt*log(2) + (Nt-1)*c),
		   "equal" = return((Nt-1)*log(2) + (Nt-1)*c + log(1+t^(Nt-1))),
		   "opposite" = return((Nt-1)*log(2) + (Nt-1)*c + log(1-t^(Nt-1))),
		   stop("No such boundaries!")
		   )
}

log.z.classical.2d = function(h, L, Nt, boundaries="periodic"){
	if(boundaries == "periodic" | boundaries == "open"){
		return(L * log.z.classical.1d(h, Nt, boundaries))
	}else{
		if(length(boundaries) == 1) boundaries = rep(boundaries, L)
		stopifnot(length(boundaries) == L)

		res = ifelse(boundaries == 1,
					 log.z.classical.1d(h, Nt, "equal"),
					 log.z.classical.1d(h, Nt, "opposite")
					 )
		return(sum(res))
	}
}

free.real.sign = function(L, Nt, h, time, ...){
	delta = time/Nt
	h_cl = -1/2*log(tanh(delta*h))
	return(exp(Re(log.z.classical.2d(h_cl, L, Nt, ...)-log.z.classical.2d(Re(h_cl), L, Nt, ...))))
}

free.real.evolution = function(L, h, time, ...){
	return((2*cos(h*Im(time)))^L)
}

norm.dos = function(L, Nt, h, time, boundaries="periodic", truncate=TRUE, quantum=TRUE, ...){
	delta = time/Nt
	h_cl = -1/2*log(tanh(delta*h))
	if(truncate) h_cl = Re(h_cl)
	if(quantum) a = 1/2*log(sinh(delta*h)*cosh(delta*h))
	else a = -h_cl
	if(boundaries == "periodic"){
		return(L*Nt*a + log.z.classical.2d(h_cl, L, Nt, boundaries))
	}else{
		return(L*log(1/2) + L*(Nt-1)*a + log.z.classical.2d(h_cl, L, Nt, boundaries))
	}
}

rescale.sff = function(sff, L, Nt, J, J2, h, time, avoid, ...){
	norm = exp(Re(norm.dos(L, Nt, h, time, ...)))
	c.norm = exp(Re(norm.dos(L, Nt, h, time, quantum=FALSE, ...)))
	offset = classical_approx(L, Nt, J, J2, h, time, avoid=avoid)
	return((sff + offset/c.norm) * norm)
}
