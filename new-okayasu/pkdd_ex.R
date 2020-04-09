require(phacm)
require(tidyverse)
require(TDA)
require(myfs)
require(rgl)
library(seephacm)


#PKDD用追加実験

#一様ノイズ付加
torus15.300sub6<-subsampleExclude(torus.collect15, nsub = 300)
save2Rdata(torus15.300sub6)

figurePlot3d(torus15.300sub6[[1]][["noizyX"]])

uninois<-cbind(runif(50, -3, 3), runif(50, -3, 3), runif(50, -1, 1))

#figurePlot3d(rbind(torus15.300sub6[[1]][["noizyX"]], uninois))
points3d(uninois, col=2)
rgl.snapshot("./uninoise.png")

noitorus15.300sub6<-torus15.300sub6
for (k in 1:100) {
  uninoise<-cbind(runif(50, -3.5, 3.5), runif(50, -3.5, 3.5), runif(50, -1, 1))
  noitorus15.300sub6[[k]][["noizyX"]]<-rbind(torus15.300sub6[[k]][["noizyX"]], uninoise)
}  

noitorus15.300insub6<-intering(noitorus15.300sub6)
points3d(noitorus15.300insub6[[1]][["noizyX"]][351:586, ], col=3)
rgl.snapshot("./ineduninoise.png")

noitorus15.300insub6.aggr<-proposedMethodOnly(X = noitorus15.300insub6, maxdim = 2, maxscale = 3, samples = 10)
save2Rdata(noitorus15.300insub6.aggr)

torus15.300sub6.aggr<-proposedMethodOnly(X = torus15.300sub6, maxdim = 2, maxscale = 3, samples = 10)
save2Rdata(torus15.300sub6.aggr)

#ガウス分布の値を足して点を散らす
gautorus15.300sub6<-torus15.300sub6
for (k in 1:100) {
  nsample<-length(torus15.300sub6[[k]][["noizyX"]][,1])
  subsam<-noitorus15.300sub6[[k]][["noizyX"]][sample(nsample, nsample*0.7), ]
  for (l in 1:length(subsam[,1])) {
    subsam[l, ]<-subsam[l, ]+rnorm(n = 3, mean = 0, sd = 0.2)
  }
  gautorus15.300sub6[[k]][["noizyX"]]<-rbind(torus15.300sub6[[k]][["noizyX"]], subsam)
  gautorus15.300sub6[[k]][["nsample"]]<-nsample+nsample*0.7
}  

figurePlot(gautorus15.300sub6[[1]][["noizyX"]][1:300,])
points3d(gautorus15.300sub6[[1]][["noizyX"]][301:510,], col=3)
rgl.snapshot("./gauss_add.png")

gautorus15.300sub6.aggr<-proposedMethodOnly(X = gautorus15.300sub6, maxdim = 2, maxscale = 3, samples = 10)
save2Rdata(gautorus15.300sub6.aggr)
