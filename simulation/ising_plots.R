source("ising_dos.R")
source("ising_plot-functions.R")

# auxiliary params
n.boot = 40
steps = c(1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8)

# color coding for plots
cols = c("violet", "blue", "green", "orange", "red")

i = 1
dir = c("short", "long")[i]
# default lattice extent
sizes = switch(i, c(4, 6, 8, 10, 12), c(4))
offset = ifelse(i == 1, 0, 2^sizes)
avoid = 2

L = 16
Nt = 4
t = 1i/min(L, 12)
n_states = ifelse(L == 6, 120, 10)

# vary t
k = offset + 1:min(L, 12)

# physical parameters
J2 = .3
h = .6
dJ = 1

J = read.table(dir(paste0("../production/L_", L, "_dJ_", dJ, "/"), pattern="*_field.txt", full.name=TRUE)[1])[,1]

t.scale = (dJ-1)*2/5 + 1
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[1])
load(files[2])

# exact results
en = read.table(dir(paste0("../production/L_", L, "_dJ_", dJ, "/"), pattern="*_ev.txt", full.name=TRUE)[1])[,1]

file =  paste0("../plots/sff_llr_L_", L, "_dJ_", dJ, "_av_", avoid)
pdf(file=paste0(file, ".pdf"), pointsize=14)
plot.dos(Im(t/t.scale*k), res.llr, steps, "N_s", "sff", ylab=TeX("$K(t)$"), log="y", ylim=c(.5,2^L)^2)
plot.exact(Im(t/t.scale*k), energies=en)
plot.classical(Im(t/t.scale*k), L, J, J2)
plot.2nd.order(Im(t/t.scale*k), L, J, J2, h)
plot.dos.spline(k, res.llr, steps, "N_s", "sff", L, Nt, J, J2, h, t/t.scale, avoid, 6, ylab=TeX("$K(t)$"), log="y", ylim=c(.5,2^L)^2)
plot.exact(Im(t/t.scale*k), energies=en)
plot.classical(Im(t/t.scale*k), L, J, J2)
plot.2nd.order(Im(t/t.scale*k), L, J, J2, h)
dev.off()

pdf(file=paste0("../plots/sff_rew_L_", L, "_dJ_", dJ, "_av_", avoid, ".pdf"), pointsize=14)
plot.dos(Im(t/t.scale*k), res.rew, steps, "N_s", "sff", ylab=TeX("$K(t)$"), log="y", ylim=c(.5,2^L)^2)
plot.exact(Im(t/t.scale*k), energies=en)
plot.classical(Im(t/t.scale*k), L, J, J2)
plot.2nd.order(Im(t/t.scale*k), L, J, J2, h)
dev.off()

pdf(file=paste0("../plots/spline_L_", L, "_dJ_", dJ, ".pdf"), pointsize=14)
n = 1
plot.spline(res.llr, n, Im(t/t.scale*k), main=paste0("L = ", L, ", dJ = ", dJ, ", steps = ", steps[n]), line.col="red", ylim=c(1e-6,1e-1)^2)
dev.off()

#pdf(file=paste0("../plots/DoS_L_", L, "_dJ_", dJ, "_av_", avoid, ".pdf"), pointsize=14)
#pos = L
#plot.probs(get.ImS(J, J2, k[pos]*t/t.scale, k[pos]*n_states), res.llr, 9, pos, ylab=TeX("$\\rho(S_I)$"), log="y", ylim=c(1e-6,1e-1))
#dev.off()

pdf(file=paste0("../plots/error_L_", L, "_dJ_", dJ, "_av_", avoid, ".pdf"), pointsize=14)
plot.transposed(Im(t/t.scale*k), res.llr, steps, "N_s", "sff", energies=en, ylab=TeX("$\\Delta K(t)$"), log="xy", ylim=c(1e-7,10)*2)
plot.stat.pow(Im(t/t.scale*k), energies=en, Nt=Nt*k)
plot.transposed.spline(k, res.llr, steps, "N_s", "sff", energies=en, L, Nt, J, J2, h, t/t.scale, avoid, 6, ylab=TeX("$\\Delta K(t)$"), log="xy", ylim=c(1e-7,10)*2)
plot.stat.pow(Im(t/t.scale*k), energies=en, Nt=Nt*k)
plot.transposed(Im(t/t.scale*k), res.rew, steps, "N_s", "sff", energies=en, ylab=TeX("$\\Delta K(t)$"), log="xy", ylim=c(1e-7,10)*2)
plot.stat.pow(Im(t/t.scale*k), energies=en, Nt=Nt*k)

#dJ = 0
#files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
#load(files[1])
#
#pdf(file=paste0("../plots/DoS_L_", L, "_dJ_", dJ, "_av_", avoid, ".pdf"), pointsize=14)
#plot(rho$ImS, rho$dos[,1], xlab=TeX("$S_I$"), ylab=TeX("$\\rho_+(S_I)$"), log="y")
#lines(rho$ImS, rho$dos[,1])
#dev.off()
#dJ = 1

avoid = 0
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[1])

plot.transposed(Im(t/t.scale*k), res.llr, steps, "N_s", "sff", energies=en, ylab=TeX("$\\Delta K(t)$"), log="xy", ylim=c(1e-7,10)*2)
plot.stat.pow(Im(t/t.scale*k), energies=en, Nt=Nt*k)
dev.off()

#pdf(file=paste0("../plots/DoS_L_", L, "_dJ_", dJ, "_av_", avoid, ".pdf"), pointsize=14)
#pos = L
#plot.probs(get.ImS(J, J2, k[pos]*t/t.scale, k[pos]*n_states), res.llr, 9, pos, ylab=TeX("$\\rho(S_I)$"), log="y", ylim=c(1e-6,1e-1))
#dev.off()

pdf(file=paste0("../plots/sff_llr_L_", L, "_dJ_", dJ, "_av_", avoid, ".pdf"), pointsize=14)
plot.dos(Im(t/t.scale*k), res.llr, steps, "N_s", "sff", ylab=TeX("$K(t)$"), log="y", ylim=c(.5,2^L)^2)
plot.exact(Im(t/t.scale*k), energies=en)
plot.classical(Im(t/t.scale*k), L, J, J2)
plot.2nd.order(Im(t/t.scale*k), L, J, J2, h)
dev.off()

n_states = 10
pdf(file=paste0("../plots/sff_llr_large_dJ_", dJ, "_av_", avoid, ".pdf"), pointsize=14)
L = 20
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[1])
plot.dos(Im(t/t.scale*k), res.llr[7], seq(20, 50, by=5), "L", "sff", scale=.5^L, ylab=TeX("$2^{-2L}K(t)$"), log="y", ylim=c(.5^10,1)^2)

L = 30
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[1])
plot.dos(Im(t/t.scale*k), res.llr[7], NULL, "L", "sff", scale=.5^L, new.plot=FALSE, col.offset=1, ylab=TeX("$2^{-2L}K(t)$"))

L = 40
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[1])
plot.dos(Im(t/t.scale*k), res.llr[1], NULL, "L", "sff", scale=.5^L, new.plot=FALSE, col.offset=2, ylab=TeX("$2^{-2L}K(t)$"))

L = 50
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[1])
plot.dos(Im(t/t.scale*k), res.llr[1], NULL, "L", "sff", scale=.5^L, new.plot=FALSE, col.offset=3, ylab=TeX("$2^{-2L}K(t)$"))
dev.off()

avoid = 2
pdf(file=paste0("../plots/sff_llr_large_dJ_", dJ, "_av_", avoid, ".pdf"), pointsize=14)
L = 20
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[1])
plot.dos(Im(t/t.scale*k), res.llr[7], seq(20, 40, by=5), "L", "sff", scale=.5^L, ylab=TeX("$2^{-2L}K(t)$"), log="y", ylim=c(.5^10,1)^2)

L = 30
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[1])
plot.dos(Im(t/t.scale*k), res.llr[7], NULL, "L", "sff", scale=.5^L, new.plot=FALSE, col.offset=1, ylab=TeX("$2^{-2L}K(t)$"))

L = 40
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[1])
plot.dos(Im(t/t.scale*k), res.llr[1], NULL, "L", "sff", scale=.5^L, new.plot=FALSE, col.offset=2, ylab=TeX("$2^{-2L}K(t)$"))
dev.off()

avoid = 0
pdf(file=paste0("../plots/sff_rew_large_dJ_", dJ, "_av_", avoid, ".pdf"), pointsize=14)
L = 20
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[2])
plot.dos(Im(t/t.scale*k), res.rew[7], seq(20, 50, by=5), "L", "sff", scale=.5^L, ylab=TeX("$2^{-2L}K(t)$"), log="y", ylim=c(.5^10,1)^2)

L = 30
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[2])
plot.dos(Im(t/t.scale*k), res.rew[7], NULL, "L", "sff", scale=.5^L, new.plot=FALSE, col.offset=1, ylab=TeX("$2^{-2L}K(t)$"))

L = 40
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[2])
plot.dos(Im(t/t.scale*k), res.rew[1], NULL, "L", "sff", scale=.5^L, new.plot=FALSE, col.offset=2, ylab=TeX("$2^{-2L}K(t)$"))

L = 50
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[2])
plot.dos(Im(t/t.scale*k), res.rew[1], NULL, "L", "sff", scale=.5^L, new.plot=FALSE, col.offset=3, ylab=TeX("$2^{-2L}K(t)$"))
dev.off()

avoid = 2
pdf(file=paste0("../plots/sff_rew_large_dJ_", dJ, "_av_", avoid, ".pdf"), pointsize=14)
L = 20
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[2])
plot.dos(Im(t/t.scale*k), res.rew[7], seq(20, 40, by=5), "L", "sff", scale=.5^L, ylab=TeX("$2^{-2L}K(t)$"), log="y", ylim=c(.5^10,1)^2)

L = 30
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[2])
plot.dos(Im(t/t.scale*k), res.rew[7], NULL, "L", "sff", scale=.5^L, new.plot=FALSE, col.offset=1, ylab=TeX("$2^{-2L}K(t)$"))

L = 40
times = ifelse(L <= 30, L, 12)
Nt = floor(4 * L / times)
t = 1i/times
k = offset + 1:ifelse(L < 40, min(times, 12), 3)
files = paste0("R.Data/", dir, "/res.", c("llr", "rew"), avoid, ".L", L, ".dJ", dJ, ".RData")
load(files[2])
plot.dos(Im(t/t.scale*k), res.rew[1], NULL, "L", "sff", scale=.5^L, new.plot=FALSE, col.offset=2, ylab=TeX("$2^{-2L}K(t)$"))
dev.off()
