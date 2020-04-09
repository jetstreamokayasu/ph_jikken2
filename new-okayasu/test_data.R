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

plot3d(torus7.50)
aspect3d("iso")

torus7.diag<-ripsDiag(torus7.50, maxdimension = 2, maxscale = 3)
plot(torus7.diag[[1]])

DiagGrid <- gridDiag(X = torus7.50, FUN = kde, h = 0.3, by = by, sublevel = FALSE, library = "Dionysus", location = TRUE, printProgress = FALSE)
Land.dim1.torus7 <- landscape(torus7.diag[[1]], dimension = 1, KK = 1, tseq.torus)
plot(tseq.torus, Land.dim1.torus7, type = "l", col=2, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 1.1))

torus7.land<-calcLandscape(torus7.diag)
plotLandscape(torus7.land)

torus8.200<-torusUnif(200, 1, 2.5)
plot3d(torus8.200)
aspect3d("iso")

torus8.diag<-ripsDiag(torus8.200, maxdimension = 2, maxscale = 3)
plot(torus8.diag[[1]])



torus8.land<-calcLandscape(torus8.diag)
plotLandscape(torus8.land)

torus.collect9 <- lapply(1:100, function(i){
  nsample <- 200
  #var <- runif(1, var.min, var.max)
  #noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, noizyX = torus, diag = 0))
})
save(torus.collect9, file = "./data/torus.collect9")

#トーラス状データのサイクル数推定
torus9.aggr<-homMethodsComp2compari3(torus.collect9, 2, 3, 10)
save(torus9.aggr, file = "./data/torus9.aggr")

torus9.esit.dim1<-cyclenumber(torus9.aggr[[1]])
torus9.esit.dim2<-cyclenumber(torus9.aggr[[2]])


torus.collect10 <- lapply(1:100, function(i){
  nsample <- 300
  #var <- runif(1, var.min, var.max)
  #noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, noizyX = torus, diag = 0))
})
save(torus.collect10, file = "./data/torus.collect10")

figurePlot(torus.collect10[[1]][["noizyX"]])

#トーラス状データのサイクル数推定
torus10.aggr<-homMethodsComp2compari3(torus.collect10, 2, 3, 10)
save(torus10.aggr, file = "./data/torus10.aggr")

torus10.esit.dim1<-cyclenumber(torus10.aggr[[1]])
torus10.esit.dim2<-cyclenumber(torus10.aggr[[2]])


#関数がうまくいくか試す用の球状データセット
sphere.collect7<- lapply(1:100, function(i){
  sphere <- sphereUnif(300, 2, 1)
  #cat("nsample=", nsample, "data", i, "\n")
  return(list(nsample = 300, noizyX = sphere, diag = 0))
})
save(sphere.collect7, file = "./data/sphere.collect7")

#球状データのサイクル数推定
sphehole7.aggr<-homMethodsComp2compari2(sphere.collect7, 2, 3, 10)
save(sphehole3.aggr, file = "./data/sphehole3.aggr")

#近傍15点の距離の平均による閾値設定
sphehole2.aggr1.land<-calcLandscape(sphehole2.aggr[["Xdiag"]][[1]])
pre.thresh<-meanVicsDestance(sphere.collect2[[1]][["noizyX"]], 15)
abline(h=pre.thresh/2, col="orange")
abline(h=pre.thresh/4, col="orange")

sphehole2.aggr2.land<-calcLandscape(sphehole2.aggr[["Xdiag"]][[2]])
pre.thresh2<-meanVicsDestance(sphere.collect2[[2]][["noizyX"]], 15)
abline(h=pre.thresh2/2, col="orange")
abline(h=pre.thresh2/4, col="orange")

sphehole2.aggr3.land<-calcLandscape(sphehole2.aggr[["Xdiag"]][[3]])
pre.thresh3<-meanVicsDestance(sphere.collect2[[3]][["noizyX"]], 15)
abline(h=pre.thresh3/2, col="orange")
abline(h=pre.thresh3/4, col="orange")

sphehole2.aggr4.land<-calcLandscape(sphehole2.aggr[["Xdiag"]][[4]])
pre.thresh4<-meanVicsDestance(sphere.collect2[[4]][["noizyX"]], 15)
abline(h=pre.thresh4/2, col="orange")
abline(h=pre.thresh4/4, col="orange")

data.collect1.land<-calcLandscape(list(diagram=data.collect[[1]][["diag"]]))
torus.thresh1<-meanVicsDestance(data.collect[[1]][["noizyX"]], 15)
abline(h=torus.thresh1/2, col="orange")
abline(h=torus.thresh1/4, col="orange")

#ノイズトーラスのサイクル数推定結果を見る
data.hole.dim1<-cyclenumber(hole.aggr2[[1]])
data.hole.dim2<-cyclenumber(hole.aggr2[[2]])

#phacmテスト##############
data.cole1.pl<-compute_pl(data.collect[[1]][["diag"]])
autoplot(data.cole1.pl)

data.cole1.smpl<-compute_smooth_pl(data.cole1.pl)
autoplot(data.cole1.smpl)

data.cole1.pl.loma<-count_local_maximal(data.cole1.pl, 0.1)
data.cole1.smpl.loma<-count_smooth_maximal(data.cole1.pl)

sphe2.pl<-compute_pl(sphehole2.aggr[["Xdiag"]][[1]])
autoplot(sphe2.pl)

sphe2.pl.loma<-count_smooth_maximal(sphe2.pl)

##########################
#推定精度が下がる密度で単位球を試す
sphere.40<-sphereUnif(41, 2, 1)
plot3d(sphere.40)
sphe.40.diag<-ripsDiag(sphere.40, 2, 3)
plot(sphe.40.diag[[1]])


sphere.collect9<- lapply(1:100, function(i){
  sphere <- sphereUnif(38, 2, 1)
  #cat("nsample=", nsample, "data", i, "\n")
  return(list(nsample = 38, noizyX = sphere, diag = 0))
})

plot3d(sphere.collect9[[1]][["noizyX"]])

sphe9.aggr<-proposedMethodOnly(sphere.collect9, 2, 3, 10, const.size = 30)
