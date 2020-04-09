# Xlim <- c(-4, 4);
# Ylim <- c(-4, 4);
# Zlim <- c(-4, 4);
# by <- 0.1;
# Xseq <- seq(Xlim[1], Xlim[2], by = by)
# Yseq <- seq(Ylim[1], Ylim[2], by = by)
# Grid <- expand.grid(Xseq, Yseq)

torus10.knn.Diag <- gridDiag(X = torus.collect10[[1]][["noizyX"]], FUN = knnDE, k = 100,
                                   lim = cbind(Xlim, Ylim, Zlim),by = by,
                                   sublevel = FALSE, library = "Dionysus",
                                   printProgress = TRUE)

plot(torus10.knn.Diag[["diagram"]])

torus10.KDE.Diag <- gridDiag(X = torus.collect10[[1]][["noizyX"]], FUN = kde, h=1,
                             lim = cbind(Xlim, Ylim, Zlim),by = by,
                             sublevel = FALSE, library = "Dionysus",
                             printProgress = TRUE)

plot(torus10.KDE.Diag[["diagram"]])


m = 100
by2=(Xlim[2]-Xlim[1])/m
Xseq = seq(Xlim[1], Xlim[2], by=by2)
Yseq = seq(Ylim[1], Ylim[2], by=by2)
mx=length(Xseq)
my=length(Yseq)
Grid = expand.grid(Xseq,Yseq)

torus10.Diag.dist = gridDiag(torus.collect10[[1]][["noizyX"]], distFct, lim=cbind(Xlim, Ylim, Zlim), by=by2,sublevel=TRUE,printProgress=TRUE)
plot(torus10.Diag.dist[["diagram"]])

torus10.Diag.dtm = gridDiag(torus.collect10[[1]][["noizyX"]], dtm, lim=cbind(Xlim, Ylim, Zlim), by=by2,sublevel=TRUE,printProgress=TRUE, m0=0.1)
plot(torus10.Diag.dtm[["diagram"]])

sphere8.diag.dist


#追加実験
Xlim <- range(data.collect[[1]][["noizyX"]][,1])
Ylim <- range(data.collect[[1]][["noizyX"]][,2])
Zlim <- range(data.collect[[1]][["noizyX"]][,3])
by <- 0.065
Xseq <- seq(Xlim[1], Xlim[2], by = by)
Yseq <- seq(Ylim[1], Ylim[2], by = by)
Zseq <- seq(Zlim[1], Zlim[2], by = by)
Grid <- expand.grid(Xseq, Yseq, Zseq)

data.collect1.DiagGrid <- gridDiag(
   X = data.collect[[1]][["noizyX"]], FUN = kde, h = 0.3, lim = cbind(Xlim, Ylim, Zlim), by = by,
   sublevel = FALSE, library = "Dionysus", location = TRUE,
   printProgress = FALSE)
plot(data.collect1.DiagGrid$diagram)

band.kde <- bootstrapBand(X = data.collect[[1]][["noizyX"]], FUN = kde, h=0.3, Grid = Grid, B = 10, parallel = F, alpha = 0.05)

plot(data.collect1.DiagGrid$diagram, band = band.kde[["width"]])
abline(band.kde[["width"]],1, col = "red", lty = 2)

bottle.boot <- boot(data=data.collect[[1]][["noizyX"]],
                    statistic=supnor.bottle,
                    R=10,
                    Diag=data.collect1.DiagGrid)
q95.bottle <- quantile(bottle.boot$t,probs = 0.95)
abline(q95.bottle,1, col = "red", lty = 2)
 
band <- bootstrapBand(X = data.collect[[1]][["noizyX"]], FUN = kde, h=0.3, Grid = Grid, B = 10,parallel = TRUE, alpha = 0.2)
plot(data.collect1.DiagGrid$diagram, band=band[["width"]])


#bottleneck bootstrap試し
bottle.boot.dim1 <- boot(data=data.collect[[1]][["noizyX"]],
                    statistic=supnor.bottle,
                    R=10,
                    Diag=data.collect1.DiagGrid[[1]][data.collect1.DiagGrid[[1]][,1]==1,])




data.collect1.kde<-kde(data.collect[[1]][["noizyX"]], Grid, h=0.3, printProgress = TRUE)



plot(data.collect1.kde[[1]])

supnor.boot <- boot(data=data.collect[[1]][["noizyX"]],
                    statistic=supnor,
                    R=10,
                    parallel="multicore",
                    ncpus=3, kde=data.collect1.ked)
q95 <- quantile(supnor.boot$t,probs = 0.95)
plot(data.collect1.DiagGrid$diagram)
abline(q95,1, col = "black", lty = 2)

data.collect1.knn.Diag <- gridDiag(X = data.collect[[1]][["noizyX"]], FUN = knnDE, k = 100,
                            lim = cbind(Xlim, Ylim, Zlim),by = by,
                            sublevel = FALSE, library = "Dionysus",
                            printProgress = TRUE)

band.knn <- bootstrapBand(X = data.collect[[1]][["noizyX"]], FUN = knnDE,k=100, Grid = Grid, B = 10, parallel = TRUE, alpha = 0.2)

plot(data.collect1.knn.Diag[[1]], band=band.knn[["width"]])

data.collect1.Diagdtm <- gridDiag(
  X = data.collect[[1]][["noizyX"]], FUN = dtm, lim = cbind(Xlim, Ylim, Zlim), by = by,
  sublevel = TRUE, m0=0.1, printProgress = TRUE)
