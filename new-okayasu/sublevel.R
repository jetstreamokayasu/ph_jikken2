require(phacm)
require(tidyverse)
require(TDA)
require(myfs)
require(rgl)
require(boot)
require(doParallel)
require(foreach)
require(plotly)
require(viridis)
require(pterrace)
require(phacm)


#sin波のsublevel setが計算できるか試す
phi<-seq(0, 4*pi, length=200)
plot(phi, sin(phi))

DiagGrid <- gridDiag(X = sin(phi), FUN = kde, h = 0.1, lim = c(-2, 2), by = 0.1,sublevel = T, library = "Dionysus", location = TRUE, 
                     printProgress = FALSE)

DiagGrid2 <- gridDiag(X = cos(phi), FUN = kde, h = 0.1, lim = c(-1, 1), by = 0.1,sublevel = T, library = "Dionysus", location = TRUE, 
                     printProgress = FALSE)

DiagGrid3 <- gridDiag(X = sin(phi)+sin(2*phi), FUN = kde, h = 0.1, lim = c(-2, 2), by = 0.1, sublevel = T, library = "Dionysus", location = TRUE, 
                      printProgress = FALSE)

DiagGrid4 <- gridDiag(X = c(cos(phi[101:200])-1, sin(phi[1:100])+sin(2*phi[1:100])), FUN = kde, h = 0.1, lim = c(-2, 1.8), by = 0.1, sublevel = T, library = "Dionysus", location = TRUE, 
                      printProgress = FALSE)


Xlim <- c(0, 10);
Ylim <- c(0, 10);
by <- 0.1;
Xseq <- seq(Xlim[1], Xlim[2], by = by)
Yseq <- seq(Ylim[1], Ylim[2], by = by)
Grid <- expand.grid(Xseq, Yseq)

KDE<-TDA::kde(sin(phi), Grid = seq(0, 13, by = by), h=0.5)
plot(seq(0, 13, by = by), KDE)

#極値を求めてみる
cos_crip<-critical_point(cos(phi))
sin_crip<-critical_point(sin(phi)+sin(2*phi))


cos_comp<-component_set(cos_crip)
sin_comp<-component_set(sin_crip)


sin_set<-sublevelset_persist(sin(phi)+sin(2*phi), crip = sin_crip)
plot(sin_set-min(sin_set), xlim = c(0, 4), ylim = c(0, 4), pch=16)
abline(0, 1)

cos_set<-sublevelset_persist(cos(phi), crip = cos_crip)
