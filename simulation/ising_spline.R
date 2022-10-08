spline.ft = function(sp){
	if(is.atomic(sp)) return(NA)

	x = seq(min(sp$x), max(sp$x), length.out = 1001)
	dx = x[2]-x[1]
	y = exp(predict(sp, x)$y + 1i*x)

	w = c(1, 4, rep(c(2, 4), (length(y)-3)/2), 1)

	return(sum(w * y)/3 * dx)
}

spline.fit = function(rho, spline.order){
	if(missing(spline.order)) spline.order = rho$spline.order
	if(is.na(spline.order)) return(NULL)

	norm = rho$renorm * (rho$topo.ratio - 1) / (rho$topo.ratio + 1)

	x = rho$ImS
	y = (rho$dos[,1] - rho$dos[,2]) / (x[2] - x[1])

	mask = which(y > 0)
	x = x[mask]
	dy = 1/sqrt(y[mask]*rho$steps)
	y = log(y[mask])

	res = smooth.spline(x, y, df=spline.order, w=1/dy^2)

	return(res)
}

spline.sff.w.err = function(ImS, dos, samples){
	n = dim(dos)[2]/2
	dx = ImS[2] - ImS[1]
	y0 = (dos[1,1:n] - dos[1,-(1:n)]) / dx
	mask = which(y0 > 0)

	x = ImS[mask]
	dy = rep(1, length(mask))
	y0 = log(y0[mask])

	i = c(1, length(dy))
	dy[i] = dy[i] / 3

	y.s = apply(samples, 2, function(y) ((y[1:n] - y[-(1:n)]) / dx)[mask])
	#now bootstrap
	n = dim(samples)[2]
	y.s = sapply(1:n, function(k){
					 kk = floor(n*runif(n))
					 apply(y.s[,kk], 1, median)
				})

	order = floor((length(x)-2)^seq(0, 1, length.out=11)) + 2

	res = lapply(order, function(df){
					 f0 = smooth.spline(x, y0, df=df, w=1/dy^2)
					 f.s = apply(y.s, 2, function(y){
									 msk = which(y > 0)
									 f = tryCatch(
												  smooth.spline(x[msk], log(y[msk]), df=df, w=1/dy[msk]^2),
												  error = function(cond){
													  message(cond)
													  message("\n")
													  return(NA)
												  })
									 return(f)
								})
					 return(list(f0=f0, f.s=f.s))
				})

	sff = sapply(res, function(sp){
					 int = spline.ft(sp$f0)
					 dist = sapply(sp$f.s, spline.ft)
					 re = sd(Re(dist), na.rm=TRUE)
					 im = sd(Im(dist), na.rm=TRUE)
					 return(c(int, re + 1i*im))
				})

	n = dim(dos)[2]/2
	dy = sqrt(apply(dos[2:3,1:n]^2 + dos[2:3,-(1:n)]^2, 2, mean)) / dx
	dy = dy[mask] / exp(y0)
	n = length(mask)

	bic = sapply(seq(res), function(i){
					 df = order[[i]]
					 sp = res[[i]]
					 chi_sq = sum(((y0 - predict(sp$f0, x)$y) / dy)^2)
					 return(n*log(chi_sq/n) + df*log(n))
				})

	return(list(ImS=ImS, res=res, sff=sff, order=order, bic=bic))
}

spline.meas = function(rho, spline.order){
	res = spline.fit(rho, spline.order)
	return(spline.ft(res))
}

plot.smooth.spline = function(sp, ImS, dos, line.col="black", ...){
	n = dim(dos)[2]/2
	dx = ImS[2] - ImS[1]
	y = (dos[1,1:n] - dos[1,-(1:n)]) / dx
	dy = sqrt(apply(dos[2:3,1:n]^2 + dos[2:3,-(1:n)]^2, 2, mean)) / dx
	plotwitherror(x=ImS, y=y, dy=dy, ...)

	x = seq(min(sp$f0$x), max(sp$f0$x), length.out = 400)
	y = exp(predict(sp$f0, x)$y)

	lines(x, y, col=line.col, lwd=3)
}
