require(TDA)
require(myfs)
require(rgl)
nsample.min <- 300
nsample.max <- 500
var.min <- 0
var.max <- 0.3
sample<-500
ndata <- 100#データセット数

maxdimension=2
maxscale=3

#ノイズを乗せたトーラス上データを100セット
data.collect <- lapply(1:ndata, function(i){
  nsample <- round(runif(1, nsample.min, nsample.max))
  var <- runif(1, var.min, var.max)
  noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, var = var, noizyX = torus + noize.torus,
              X = torus, diag = 0))
})

#関数がうまくいくかテスト用のノイズを乗せたトーラス上データを3セット
test.collect2 <- lapply(1:3, function(i){
  nsample <- round(runif(1, nsample.min, nsample.max))
  var <- runif(1, var.min, var.max)
  noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, var = var, noizyX = torus + noize.torus,
              X = torus, diag = 0))
})

#テスト用データのサイクル数推定
test2.aggr<-homMethodsComp2compari(test.collect2, 2, 3, 10)
plot(test2.aggr[["Xdiag"]][[1]]$diagram)

#save(data.collect, file = "./data/data-collect")

#生成したトーラスデータのPH図を算出
for (i in 1:ndata) {
  data.collect[[i]]$diag <- ripsDiag(data.collect[[i]]$noizyX, maxdimension = maxdimension, maxscale = maxscale)[[1]]
}

#save(data.collect, file = "./data/data-collect")

#ノイズを乗せたトーラスのサイクル数推定
hole.aggr<-homMethodsComp2compari(data.collect, 2, 3, 10)
hole.aggr2<-homMethodsComp2compari(data.collect, 2, 3, 10)
#save(hole.aggr, file = "./data/hole-aggr")
#save(hole.aggr2, file = "./data/hole-aggr2")

plot(data.collect[[1]]$diag)

plot3d(data.collect[[2]]$noizyX, add = T)

plot3d(data.collect[[2]]$noizyX)
aspect3d("iso")

# test.collect <- lapply(1:3, function(i){
#   nsample <- round(runif(1, nsample.min, nsample.max))
#   var <- runif(1, var.min, var.max)
#   noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
#   torus <- torusUnif(nsample, 1, 2.5)
#   return(list(nsample = nsample, var = var, noizyX = torus + noize.torus,
#               X = torus, diag = 0))
# })

#test.aggr<-homMethodsComp2compari(test.collect, 2, 3, 10)
test.aggr2<-homMethodsComp2compari(test.collect, 2, 3, 10)

#それぞれのデータセットの推定サイクル数をプロット
#横軸データ点数、縦軸ノイズの分散
A<-plot.result(3, 1, 2,test.collect, test.aggr)
A2<-plot.result(3, 2, 1,test.collect, test.aggr)
B<-plot.result(3, 2, 1,test.collect, test.aggr2)

#サイクル数推定結果を描写
hole.dim1<-plot.result(ndata, 1, 2, data.collect, hole.aggr)
hole.dim2<-plot.result(ndata, 2, 1, data.collect, hole.aggr)

hol2.dim1<-plot.result(ndata, 1, 2, data.collect, hole.aggr2)
hol2.dim2<-plot.result(ndata, 2, 1, data.collect, hole.aggr2)

#間違ったサイクル数を推定されたデータセットのPH図を描写
plotWrongDiag(data.collect, hole.dim1, "proposed", 2)
plotWrongDiag(test.collect, B, "proposed", 1)

#球状データセットを作成
sphere.collect <- lapply(1:ndata, function(i){
  torus <- sphereUnif(sample, 2, 1)
  return(list(nsample = sample, noizyX = torus, diag = 0))
})

#球状データのサイクル数推定
sphehole.aggr<-homMethodsComp2compari(sphere.collect, 2, 3, 10)

plot3d(sphere.collect[[1]]$noizyX)

sphere.collect[[1]]$diag<- ripsDiag(sphere.collect[[1]]$noizyX, maxdimension = maxdimension, maxscale = maxscale)[[1]]

plot(sphere.collect[[1]]$diag)

#球状データセットを作成
sphere.collect2<- lapply(1:ndata, function(i){
  nsample <- round(runif(1, nsample.min, nsample.max))
  sphere <- sphereUnif(nsample, 2, 1)
  #cat("nsample=", nsample, "data", i, "\n")
  return(list(nsample = nsample, noizyX = sphere, diag = 0))
})

plot3d(sphere.collect2[[1]]$noizyX)

sphere.collect2[[1]]$diag<- ripsDiag(sphere.collect2[[1]]$noizyX, maxdimension = maxdimension, maxscale = maxscale)[[1]]

save(sphere.collect2, file = "./data/sphere.collect2")

plot(sphere.collect2[[1]]$diag)

sphe <- sphereUnif(100, 2, 0.5)
plot3d(sphe)

#球状データのサイクル数推
sphehole2.aggr<-homMethodsComp2compari(sphere.collect2, 2, 3, 10)
save(sphehole2.aggr, file = "./data/sphehole2.aggr")

for (i in 1:ndata) {
  sphere.collect2[[i]]$diag <- ripsDiag(sphere.collect2[[i]]$noizyX, maxdimension = maxdimension, maxscale = maxscale)[[1]]
  cat("sphere data", i, "calculating\n")
}

#サイクルの有無判定
sphehole2.hole<-holeexistdiscri(sphere.collect2)

hole.exist.data<-which(sphehole2.hole[,2]==1)

#サイクル数をn個と判定されたデータセットがそれぞれいくつあるかカウント
sphere2.dim1.cycle<-cyclenumber(sphehole2.aggr[[1]])
sphere2.dim2.cycle<-cyclenumber(sphehole2.aggr[[2]])

for(l in 1:length(hole.exist.data)){
  
  cat("data", hole.exist.data[l], "=", sphere.collect2[[hole.exist.data[l]]]$nsample, "\n")
  
}

#トーラス上データセットを作成
#2次のパーシステンスに関する閾値を変えた後での推定制度を調べる用
torus.collect <- lapply(1:ndata, function(i){
  nsample <- sample
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, noizyX = torus, diag = 0))
})

save(torus.collect, file = "./data/torus.collect")

plot3d(torus.collect[[1]]$noizyX)
aspect3d("iso")

#トーラス状データのサイクル数推定
torus.aggr<-homMethodsComp2compari(torus.collect, 2, 3, 10)
save(torus.aggr, file = "./data/torus.aggr")

#サイクルの有無判定
#及びサイクル数がn個あると判されたデータセットがそれぞれいくつあるか調べる
torus.holes<-holexistdetermine(torus.aggr[["Xdiag"]])
esti.dim1.torus<-cyclenumber(torus.aggr[[1]])
esti.dim2.torus<-cyclenumber(torus.aggr[[2]])

sphehole2.hole2<-holexistdetermine(sphehole2.aggr[["Xdiag"]])
esti.dim1.sphere<-cyclenumber(sphehole2.aggr[[1]])
esti.dim2.sphere<-cyclenumber(sphehole2.aggr[[2]])

esti2.dim1.sphere<-cyclenumber2(sphehole2.aggr[[1]], sphehole2.hole2[["holeexist"]], 1)
esti2.dim2.sphere<-cyclenumber2(sphehole2.aggr[[2]], sphehole2.hole2[["holeexist"]], 2)

#サイクル・ノイズ判定の閾値を穴の次元ごとに求めた平均にする試し
test.ret<-calcDiagCentroid.mk1(torus.aggr[["Xdiag"]][[1]]) 
test2.ret<-calcDiagCentroid.mk1(sphehole2.aggr[["Xdiag"]][[1]]) 
test3.ret<-calcDiagCentroid(data.collect[[1]]$diag)

#関数がうまくいくか試す用の球状データセット
sphere.collect3<- lapply(1:3, function(i){
  sphere <- sphereUnif(400, 2, 1)
  #cat("nsample=", nsample, "data", i, "\n")
  return(list(nsample = 400, noizyX = sphere, diag = 0))
})
save(sphere.collect3, file = "./data/sphere.collect3")

#球状データのサイクル数推定
sphehole3.aggr<-homMethodsComp2compari(sphere.collect3, 2, 3, 10)
save(sphehole3.aggr, file = "./data/sphehole3.aggr")

#サイクルの有無判定
#及びサイクル数がn個あると判されたデータセットがそれぞれいくつあるか調べる
sphehole3.hole<-holexistdetermine(sphehole3.aggr[["Xdiag"]])
esti.dim1.sphehole3<-cyclenumber2(sphehole3.aggr[[1]], sphehole3.hole[["holeexist"]], 1)
esti.dim2.sphehole3<-cyclenumber2(sphehole3.aggr[[2]], sphehole3.hole[["holeexist"]], 2)

#穴の分布が指数分布だと仮定した閾値を求める
sphere3.thresh<-threshDetermine(sphehole3.aggr[["Xdiag"]][[1]],3)
#指数分布による閾値によるサイクルの有無判定
sphehole3.hole2<-holexistdetermine2(sphehole3.aggr[["Xdiag"]])

#球状データの1次パーシステントランドスケープを作成
maxscale <- 3
tseq <- seq(0, maxscale, length = 1000) #domain
Land <- landscape(sphehole3.aggr[["Xdiag"]][[1]][[1]], dimension = 1, KK = 1, tseq)
plot(tseq, Land, type = "l", col=2, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 0.4), main ="1-degree landscape of sphere")
par(new=T)
#パーシステントランドスケープに指数分布に基づいて設定した閾値を描写
plot(tseq, rep(0.3873508/2, 1000), type = "l", xlab = "",ylab = "", ylim=c(0, 0.4))
#Sil <- silhouette(DiagRips[["diagram"]], p = 1, dimension = 1, tseq)
#球状データの2次パーシステントランドスケープを作成
Land.dim2 <- landscape(sphehole3.aggr[["Xdiag"]][[1]][[1]], dimension = 2, KK = 1, tseq)
plot(tseq, Land.dim2, type = "l", col=3, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 0.6), main ="2-degree landscape of sphere")
par(new=T)
#パーシステントランドスケープに指数分布に基づいて設定した閾値を描写
plot(tseq, rep(0.3873508/4, 1000), type = "l", xlab = "",ylab = "", ylim=c(0, 0.6))

#球状データの各次元のパーシステンスを求める
sphere3.per<-calcPersistence(sphehole3.aggr[["Xdiag"]][[1]][[1]])
sphere3.per.dim1<-calcper(sphehole3.aggr[["Xdiag"]][[1]][[1]], 1)
sphere3.per.dim2<-calcper(sphehole3.aggr[["Xdiag"]][[1]][[1]], 2)
sphere3.per.all<-c(sphere3.per.dim1, sphere3.per.dim2)
log.scale<-scale(log(sphere3.per.all))

#球状データの1次の穴のパーシステンスをに基づくヒストグラムを作成
hist(sphere3.per.dim1, breaks=seq(0,3,0.05), main="Histogram", xlab="range",  xlim=c(0,3), ylim=c(0, 60), col="#993435")
hist(sphere3.per.dim1)
par(new=T)
#ヒストグラムに指数分布に基づいて設定した閾値を描写
plot(rep(0.3873508/3, 100), seq(0, 80, length= 100), type = "l", axes=F, ann = F, xlim = c(0, 3))

#球状データの全次元の穴のパーシステンスをに基づくヒストグラムを作成
hist(exp(c(sphere3.per.dim1, sphere3.per.dim2)), breaks=seq(1,3,0.05), main="Histogram", xlab="range",  xlim=c(1,3), ylim=c(0, 60), col="#993435")
par(new=T)
#ヒストグラムに指数分布に基づいて設定した閾値を描写
plot(rep(exp(0.3873508), 100), seq(0, 80, length= 100), type = "l", axes=F, ann = F, xlim = c(1, 3),  ylim=c(0, 60))


hist(log(c(sphere3.per.dim1, sphere3.per.dim2)), breaks=seq(-10, 5,0.5), main="Histogram", xlab="range", col="#993435")
hist(log.scale, breaks=seq(-4, 4, 0.1), main="Histogram", xlab="range", col="#993435")


#関数がうまくいくか試す用のトーラス状データを3セット作成
test.collect3 <- lapply(1:3, function(i){
  nsample <- 500
  #var <- runif(1, var.min, var.max)
  #noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, noizyX = torus, diag = 0))
})
save(test.collect3, file = "./data/test.collect3")

#トーラス状データのサイクル数推定
test3.aggr<-homMethodsComp2compari(test.collect3, 2, 3, 10)
save(test3.aggr, file = "./data/test3.aggr")

#サイクルの有無判定
#及びサイクル数がn個あると判されたデータセットがそれぞれいくつあるか調べる
test3.hole<-holexistdetermine(test3.aggr[["Xdiag"]])
esti2.dim1.testorus<-cyclenumber2(test3.aggr[[1]], test3.hole[["holeexist"]], 1)
esti2.dim2.testorus<-cyclenumber2(test3.aggr[[2]], test3.hole[["holeexist"]], 2)

#トーラス状データの各次元のパーシステンスを求める
testorus3.per.dim1<-calcper(test3.aggr[["Xdiag"]][[1]][[1]], 1)
testorus3.per.dim2<-calcper(test3.aggr[["Xdiag"]][[1]][[1]], 2)
#トーラス状データの1次の穴のパーシステンスをに基づくヒストグラムを作成
hist(testorus3.per.dim1, breaks=seq(0,3,0.05), main="Histogram of torus hole", xlab="persistent range",  xlim=c(0,3), col="#993435")
testorus3.allper<-c(testorus3.per.dim1, testorus3.per.dim2*2)
#球状データの全次元の穴のパーシステンスをに基づくヒストグラムを作成
hist(testorus3.allper, breaks=seq(0,3,0.05), main="Histogram of torus hole", xlab="persistent range",  xlim=c(0,3), ylim=c(0, 80), col="#993435")
par(new=T)
plot(rep(0.2051808, 100), seq(0, 80, length= 100), type = "l", axes=F, ann = F, xlim = c(0, 3))
curve(dexp(x, length(testorus3.allper)/3), from = 0, to=3, ylim=c(0, 80), ann=F)

hist(log(testorus3.allper*10), breaks=seq(0,3,0.05), main="Histogram of torus hole", xlab="persistent range",  xlim=c(0,3), ylim=c(0, 80), col="#993435")

hist(exp(testorus3.allper))

#トーラス状データの1次パーシステントランドスケープ
tseq.torus <- seq(0, maxscale, length = 1000) #domain
Land.dim1.torus <- landscape(test3.aggr[["Xdiag"]][[1]][[1]], dimension = 1, KK = 1, tseq.torus)
plot(tseq.torus, Land.dim1.torus, type = "l", col=2, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 1.1))
par(new=T)
plot(tseq, rep(0.2679733/2, 1000), type = "l", xlab = "",ylab = "", ylim=c(0, 1.1))
#Sil <- silhouette(DiagRips[["diagram"]], p = 1, dimension = 1, tseq)

#トーラス状データの2次パーシステントランドスケープ
Land.dim2.torus <- landscape(test3.aggr[["Xdiag"]][[1]][[1]], dimension = 2, KK = 1, tseq.torus)
plot(tseq.torus, Land.dim2.torus, type = "l", col=3, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 0.3))
par(new=T)
plot(tseq, rep(0.2679733/4, 1000), type = "l", xlab = "",ylab = "", ylim=c(0, 0.3))

#トーラス状データのサイクルの有無判定
test3.hole2<-holexistdetermine2(test3.aggr[["Xdiag"]])

#単位円内に一様分布するデータ点を作成
disc<-uniDisMake(200, 1)
plot(disc)
disc.diag<-ripsDiag(disc, maxdimension = 2, maxscale = 3)
plot(disc.diag[[1]])
#単位円内に一様分布するデータ点のサイクル有無判定
disc.hole<-holexistdetermine2(disc.diag[[1]])

#ノイズを乗せたトーラス状データのサイクルの有無判定
datacollect.hole<-holexistdetermine2(data.collect[[1]][["diag"]])

#トーラス状データの1次パーシステントランドスケープ
Land.dim1.torus3 <- landscape(test3.aggr[["Xdiag"]][[3]][[1]], dimension = 1, KK = 1, tseq.torus)
plot(tseq.torus, Land.dim1.torus3, type = "l", col=2, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", main ="1-degree landscape of torus", ylim=c(0, 1.1))
par(new=T)
plot(tseq, rep(0.2051808/2, 1000), type = "l", xlab = "",ylab = "", ylim=c(0, 1.1))

#トーラス状データの2次パーシステントランドスケープ
Land.dim2.torus3 <- landscape(test3.aggr[["Xdiag"]][[3]][[1]], dimension = 2, KK = 1, tseq.torus)
plot(tseq.torus, Land.dim2.torus3, type = "l", col=3, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 0.4), main ="2-degree landscape of torus")
par(new=T)
plot(tseq, rep(0.2051808/4, 1000), type = "l", xlab = "",ylab = "", ylim=c(0, 0.4))

#指数分布に基づいた閾値設定によるサイクル数推定
torus.collect2 <- lapply(1:ndata, function(i){
  nsample <- sample
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, noizyX = torus, diag = 0))
})

save(torus.collect2, file = "./data/torus.collect2")

plot3d(torus.collect2[[1]]$noizyX)
aspect3d("iso")

torus.aggr2<-homMethodsComp2compari(torus.collect2, 2, 3, 10)
save(torus.aggr2, file = "./data/torus.aggr2")

torus2.hole<-holexistdetermine2(torus.aggr2[["Xdiag"]])
esti.dim1.torus2<-cyclenumber2(torus.aggr2[[1]], torus2.hole[["holeexist"]], 1)
esti.dim2.torus2<-cyclenumber2(torus.aggr2[[2]], torus2.hole[["holeexist"]], 2)


sphere4.500<-sphereUnif(500, 2, 1)
sphere4.diag<-ripsDiag(sphere4.500, maxdimension = 2, maxscale = 3)

sphere5.30<-sphereUnif(30, 2, 1)
plot3d(sphere5.30)
sphere5.diag<-ripsDiag(sphere5.30, maxdimension = 2, maxscale = 3)
plot(sphere5.diag[[1]])

sphere5.hole<-holexistdetermine2(sphere5.diag)

maxscale <- 3
tseq <- seq(0, maxscale, length = 1000) #domain
Land.sphere4 <- landscape(sphere5.diag[[1]], dimension = 1, KK = 1, tseq)
plot(tseq, Land.sphere4, type = "l", col=2, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 0.4), main ="1-degree landscape of sphere")
par(new=T)
plot(tseq, rep(0.4362689/2, 1000), type = "l", xlab = "",ylab = "", ylim=c(0, 0.4))

Land.sphere4.dim2 <- landscape(sphere5.diag[[1]], dimension = 2, KK = 1, tseq)
plot(tseq, Land.sphere4.dim2, type = "l", col=3, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 0.15), main ="2-degree landscape of sphere")
par(new=T)
plot(tseq, rep(0.4362689/4, 1000), type = "l", xlab = "",ylab = "", ylim=c(0, 0.15))

torus3.120<-torusUnif(120, 1, 2.5)
plot3d(torus3.120)
aspect3d("iso")
torus3.diag<-ripsDiag(torus3.120, maxdimension = 2, maxscale = 3)
plot(torus3.diag[[1]])

Land.dim1.torus3 <- landscape(torus3.diag[[1]], dimension = 1, KK = 1, tseq.torus)
plot(tseq.torus, Land.dim1.torus3, type = "l", col=2, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 1.1))

sphere.collect6<- lapply(1:ndata, function(i){
  sphere <- sphereUnif(30, 2, 1)
  #cat("nsample=", nsample, "data", i, "\n")
  return(list(nsample = 30, noizyX = sphere, diag = 0))
})
save(sphere.collect6, file = "./data/sphere.collect6")

sphere6.aggr<-homMethodsComp2compari3(sphere.collect6, 2, 3, 10)
save(sphere6.aggr, file = "./data/sphere6.aggr")
