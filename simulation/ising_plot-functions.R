require("hadron")
require(latex2exp)

plot.dos = function(x, rho, var, var.name, obs, scale=1, new.plot=TRUE, col.offset=0, ...){
	lapply(seq(cols), function(i){
			   j = 2*i-1
			   n = length(x)
			   if(length(rho) >= j)
				   plotwitherror(x=x, y=(rho[[j]][[obs]][1,1:n]*scale)^2, dy=2*rho[[j]][[obs]][1,1:n]*rho[[j]][[obs]][3,1:n]*scale^2, mdy=2*rho[[j]][[obs]][1,1:n]*rho[[j]][[obs]][2,1:n]*scale^2, rep=(!new.plot | i>1), col=cols[i+col.offset], pch=i+col.offset, xlab="t", xlim=c(min(x), max(x)), ...)
	})
	if(!is.null(var)){
		ii = 1:min(length(cols), (length(var)+1)/2)
		legend("left", col=cols[ii], pch=ii, legend=sapply(ii*2-1, function(i) TeX(sprintf("$%s=%g$", var.name, var[[i]]))))
	}
}

plot.dos.spline = function(k, rho, var, var.name, obs, L, Nt, J, J2, h, t, avoid, order, ...){
	ii = seq(k)
	y = sapply(rho, function(r){
				   sapply(ii, function(i){
							  order = which.min(r$spline[[i]]$bic)-1
							  abs(rescale.sff(r$spline[[i]]$sff[1,order], L, k[i]*Nt, J, J2, h, k[i]*t, avoid))
							})
				})
	dy = sapply(rho, function(r){
				   sapply(ii, function(i){
							  order = which.min(r$spline[[i]]$bic)-1
							  abs(rescale.sff(r$spline[[i]]$sff[2,order], L, k[i]*Nt, J, J2, h, k[i]*t, avoid=0))
							})
				})
	lapply(seq(cols), function(i){
			   j = 2*i-1
			   if(length(rho) >= j)
				   plotwitherror(x=Im(k*t), y=y[,j]^2, dy=2*y[,j]*dy[,j], rep=i>1, col=cols[i], pch=i, xlab="t", xlim=c(min(Im(k*t)), max(Im(k*t))), ...)
	})
	if(!is.null(var)){
		ii = 1:min(length(cols), (length(var)+1)/2)
		legend("left", col=cols[ii], pch=ii, legend=sapply(ii*2-1, function(i) TeX(sprintf("$%s=%g$", var.name, var[[i]]))))
	}
}

plot.probs = function(x, rho, var, obs, ...){
	j = var
	lapply(0:1, function(i){
			   range = dim(rho[[j]]$dos[[obs]])[2]/2
			   range = i*range + (1:range)
			   plotwitherror(x=x, y=rho[[j]]$dos[[obs]][1,range], dy=rho[[j]]$dos[[obs]][3,range], mdy=rho[[j]]$dos[[obs]][2,range], rep=i>0, col=cols[2*i+1], pch=i+1, xlab=TeX("$S_I$"), ...)
	})
	legend("bottom", col=cols[c(1,3)], pch=1:2, legend=c("positive", "negative"))
}

plot.diff = function(x, rho, var, obs, ...){
	j = var
	range = dim(rho[[j]][[obs]])[2]/2
	even = (1:range)
	odd = range + (1:range)
	for (sign in c(1, -1)){
		plotwitherror(x=x, y=sign*(rho[[j]][[obs]][1,even]-rho[[j]][[obs]][1,odd]), dy=sqrt(rho[[j]][[obs]][3,even]^2+rho[[j]][[obs]][3,odd]^2), rep=sign<0, col=cols[3*(sign+1)/2+2], pch=(sign+1)/2+1, xlab="Im(S)", ...)
	}
	legend("bottom", col=cols[c(5,2)], pch=2:1, legend=c("even-odd", "odd-even"))
}

plot.hermit = function(rho, n, t, main, ...){
	lapply(seq(rho[[n]]$hermit), function(k){
			   r = rho[[n]]$hermit[[k]]
			   lapply(seq(r$res), function(i){
						  #summary(r$res[[i]])
						  if(!is.atomic(r$res[[i]])){
							  plot(r$res[[i]], xlab="Im(S)", ylab=TeX("$\\rho(S_I)$"), log="y", main=paste0(main, ", t = ", t[k], ", order = ", i-1), ...)
						  }
				})
	})
}

plot.spline = function(rho, n, t, main, ...){
	if(is.null(rho[[n]]$spline)) return(NULL)
	lapply(seq(rho[[n]]$spline), function(k){
			   r = rho[[n]]$spline[[k]]
			   lapply(seq(r$res), function(i){
						  if(!is.atomic(r$res[[i]])){
							  #main = paste0(main, ", t = ", t[k], ", dof = ", r$order[[i]], ", BIC = ", r$bic[[i]])
							  main = paste0("dof = ", r$order[[i]], ", BIC = ", r$bic[[i]])
							  plot.smooth.spline(r$res[[i]], r$ImS, rho[[n]]$dos[[k]], xlab=TeX("$S_I$"), ylab=TeX("$\\rho(S_I)$"), log="y", main=main, ...)
						  }
				})
	})
}

plot.exact = function(t, energies, rep=TRUE, ...){
	field = seq(min(t)-.5, max(t), length=400)
	if(!rep) plot(field, sapply(field, function(x) abs(sum(exp(1i*x*energies)))), col="white", ...)
	lines(field, sapply(field, function(x) abs(sum(exp(1i*x*energies))))^2)
	legend("bottomleft", lty=1, legend="exact")
}

plot.classical = function(t, L, J, J2, ...){
	energies = brute.force(L, J, J2)
	field = seq(min(t)-.5, max(t), length=400)
	lines(field, sapply(field, function(x) abs(sum(exp(1i*x*energies))))^2, col="brown", lty=2)
	legend("bottom", col="brown", lty=2, legend="no flips")
}

plot.2nd.order = function(t, L, J, J2, h, ...){
	energies = brute.force(L, J, J2, split=TRUE)
	en = apply(energies, 2, sum)
	field = seq(min(t)-.5, max(t), length=400)
	one.flip = sapply(field, function(x){
						  phases = apply(exp(-2i*x*energies)*sin(2*x*energies)/(2*x*energies), 2, sum)
						  return(abs(sum(exp(1i*x*en) * (1 - .5*(x*h)^2*phases))))
	})
	#no.overlap = sapply(field, function(x){
	#						phases = apply(exp(-2i*x*energies)*sin(2*x*energies)/(2*x*energies), 2, sum)
	#						return(abs(sum(exp(1i*x*en) / (1 + (1-cos(x*h))*phases))))
	#})
	full.cos = sapply(field, function(x){
							return(abs(sum(exp(1i*x*en) * cos(x*h)^L)))
	})

	lines(field, one.flip^2, col="darkgrey", lty=4)
	#lines(field, no.overlap^2, col="darkblue", lty=5)
	lines(field, full.cos^2, col="darkblue", lty=5)
	legend("bottomright", col=c("darkgrey", "darkblue"), lty=4:5, legend=c("1 flip pair", TeX("$O(t^2)$ Trotter")))
}



plot.stat.pow = function(t, energies, Nt, off=0, ...){
	ii = 1:min(length(cols), length(t))
	iii = floor(ii*length(t)/length(ii))
	lapply(ii, function(k){
			   i = iii[k]
			   lines(c(min(steps), max(steps)), rep(abs(sum(exp(1i*t[i]*energies))) / Re(exp(norm.dos(L, Nt[i], h, 1i*t[i]))), 2), col=cols[k+off])
			   #lines(c(min(steps), max(steps)), rep(free.real.sign(L, Nt, h, 1i*t[i]), 2), col=cols[i])
	})
	#legend("bottomleft", col=cols, lty=1, legend="stat. pow.")
}

plot.transposed = function(x, rho, var, var.name, obs, energies, ...){
	ii = 1:min(length(cols), length(x))
	iii = floor(ii*length(x)/length(ii))
	y0 = sapply(x, function(t) abs(sum(exp(1i*t*energies))))
	y = sapply(rho, function(r) r[[obs]][1,])
	dy = sapply(rho, function(r) r[[obs]][3,])
	mdy = sapply(rho, function(r) r[[obs]][2,])
	lapply(ii, function(k){
			   i = iii[k]
			   plotwitherror(x=var, y=abs((y[i,]/y0[i])^2-1), dy=2*y[i,]/y0[i]*ifelse(y[i,]-y0[i]>0, dy[i,], mdy[i,])/y0[i], mdy=ifelse(y[i,]-y0[i]>0, mdy[i,], dy[i,])/y0[i], rep=k>1, col=cols[k], pch=k, xlab=TeX(sprintf("$%s$", var.name)), ...)
	})
	legend("bottomleft", col=cols, pch=seq(cols), legend=sapply(x[iii], function(t) sprintf("t=%g", t)))
}

plot.transposed.hermit = function(k, rho, var, var.name, obs, energies, L, Nt, J, J2, h, t, avoid, order, ...){
	ii = 1:min(length(cols), length(k))
	iii = floor(ii*length(k)/length(ii))
	y0 = sapply(iii, function(i) abs(sum(exp(k[i]*t*energies))))
	y = sapply(rho, function(r){
				   sapply(iii, function(i)
						  abs(rescale.sff(r$hermit[[i]]$sff[[order]][1], L, k[i]*Nt, J, J2, h, k[i]*t, avoid))
						  )
				})
	dy = sapply(rho, function(r){
				   sapply(iii, function(i)
						  abs(rescale.sff(r$hermit[[i]]$sff[[order]][2], L, k[i]*Nt, J, J2, h, k[i]*t, avoid=0))
						  )
				})
	lapply(ii, function(i){
			   plotwitherror(x=var, y=abs(y[i,]-y0[i])/y0[i], dy=dy[i,]/y0[i], rep=i>1, col=cols[i], pch=i, xlab=TeX(sprintf("$%s$", var.name)), ...)
	})
	legend("bottomleft", col=cols, pch=seq(cols), legend=sapply(ii, function(i) TeX(sprintf("$t=%g$", Im(k[iii[i]]*t)))))
}

plot.transposed.spline = function(k, rho, var, var.name, obs, energies, L, Nt, J, J2, h, t, avoid, order, ...){
	ii = 1:min(length(cols), length(k))
	iii = floor(ii*length(k)/length(ii))
	y0 = sapply(iii, function(i) abs(sum(exp(k[i]*t*energies))))
	y = sapply(rho, function(r){
				   sapply(iii, function(i){
							  order = which.min(r$spline[[i]]$bic)-1
							  abs(rescale.sff(r$spline[[i]]$sff[1,order], L, k[i]*Nt, J, J2, h, k[i]*t, avoid))
						  })
				})
	dy = sapply(rho, function(r){
				   sapply(iii, function(i){
							  order = which.min(r$spline[[i]]$bic)-1
							  abs(rescale.sff(r$spline[[i]]$sff[2,order], L, k[i]*Nt, J, J2, h, k[i]*t, avoid=0))
						  })
				})
	lapply(ii, function(i){
			   plotwitherror(x=var, y=abs((y[i,]/y0[i])^2-1), dy=2*y[i,]/y0[i]*dy[i,]/y0[i], rep=i>1, col=cols[i], pch=i, xlab=TeX(sprintf("$%s$", var.name)), ...)
	})
	legend("bottomleft", col=cols, pch=seq(cols), legend=sapply(ii, function(i) TeX(sprintf("$t=%g$", Im(k[iii[i]]*t)))))
}

plot.compare = function(x, rho.llr, rho.rew, var, var.name, obs, energies, ...){
	i = 4
	y0 = sapply(x[seq(cols)], function(t) abs(sum(exp(1i*t*energies))))
	for(j in 1:2){
		rho = list(rho.llr, rho.rew)[[j]]
		y = sapply(rho, function(r) r[[obs]][1,])
		dy = sapply(rho, function(r) r[[obs]][3,])
		mdy = sapply(rho, function(r) r[[obs]][2,])
		print(y)
		plotwitherror(x=var[seq(y[i,])], y=abs(y[i,]-y0[i])/y0[i], dy=ifelse(y[i,]-y0[i]>0, dy[i,], mdy[i,])/y0[i], mdy=ifelse(y[i,]-y0[i]>0, mdy[i,], dy[i,])/y0[i], rep=j>1, col=cols[j], pch=j, xlab=var.name, ...)
	}
	legend("bottomleft", col=cols, pch=seq(cols), legend=c("LLR", "reweighting"))
}
