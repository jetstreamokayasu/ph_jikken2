source("http://aoki2.si.gunma-u.ac.jp/R/src/weibull_par.R", encoding="euc-jp")
require(ks)
require(MaskJointDensity)

#ワイブル分布に基づいた閾値設定の試し
testorus.wei<-weibull.par(testorus3.allper, epsilon=1e-7)
curve(dweibull(x, shape = testorus.wei[1], scale = testorus.wei[2]), from = 0, to=3)
plot(x=seq(0.1, 3.1, 0.1), dweibull(x, shape = testorus.wei[1], scale = testorus.wei[2]), type='l')
x<-seq(0.1, 3.1, 0.1)

curve(pweibull(x, shape = testorus.wei[1], scale = testorus.wei[2]), from = 0, to=3)
thresh.candi<-qweibull(0.95, shape = testorus.wei[1], scale = testorus.wei[2])
hist(testorus3.allper, breaks=seq(0,3,0.05), main="Histogram of torus hole", xlab="persistent range",  xlim=c(0,3), ylim=c(0, 80), col="#993435")
par(new=T)
plot(rep(thresh.candi, 100), seq(0, 80, length= 100), type = "l", axes=F, ann = F, xlim = c(0, 3))


sphere3.wei<-weibull.par(sphere3.per.all, epsilon=1e-7)
shethresh.candi<-qweibull(0.95, shape = sphere3.wei[1], scale = sphere3.wei[2])
hist(sphere3.per.dim1, breaks=seq(0,3,0.05), main="Histogram", xlab="range",  xlim=c(0,3), ylim=c(0, 60), col="#993435")
par(new=T)
plot(rep(shethresh.candi, 100), seq(0, 80, length= 100), type = "l", axes=F, ann = F, xlim = c(0, 3))

sphere3.per.dim0<-calcper(sphehole3.aggr[["Xdiag"]][[1]][[1]], 0)
hist(sphere3.per.dim0, breaks=seq(0,3,0.05), main="Histogram", xlab="range",  xlim=c(0,3), ylim=c(0, 160), col="#993435")

testorus3.dim0per<-calcper(test3.aggr[["Xdiag"]][[1]][[1]], 0)
hist(testorus3.dim0per, breaks=seq(0,3,0.05), main="Histogram of torus hole", xlab="persistent range",  xlim=c(0,3), ylim=c(0, 80), col="#993435")




sphere5.per.dim1<-calcper(sphere5.diag[[1]], 1)
sphere5.per.dim2<-calcper(sphere5.diag[[1]], 2)
sphere5.per.all<-c(sphere5.per.dim1, sphere5.per.dim2)
sphere5.wei<-weibull.par(testorus3.allper, epsilon=1e-7)
sphe5thresh.candi<-qweibull(0.95, shape = sphere5.wei[1], scale = sphere5.wei[2])
hist(sphere5.per.all, breaks=seq(0,3,0.05), main="Histogram of torus hole", xlab="persistent range",  ylim=c(0,5), col="#993435")
par(new=T)
plot(rep(sphe5thresh.candi, 100), seq(0, 80, length= 100), type = "l", axes=F, ann = F, xlim = c(0, 3))
curve(dweibull(x, shape = sphere5.wei[1], scale = sphere5.wei[2]), from = 0, to=3, ylim=c(0,5))

#KS検定
test.result<-ks.test(testorus3.allper, "pweibull", shape = testorus.wei[1], scale = testorus.wei[2], exact = T)
tst.result.sphe3<-ks.test(sphere3.per.all, "pweibull", shape = sphere3.wei[1], scale = sphere3.wei[2], exact = T)

#KDEによる閾値設定試し
sphe2.per<-calcPerdim1dim2(sphere.collect2[[35]][["diag"]])
hist(sphe2.per, breaks=seq(0,3,0.05), main="Histogram of sphere hole", xlab="persistent range", ylim=c(0,100), col="#993435")
sphe2.perdim2<-calcper(sphere.collect2[[35]][["diag"]], 2)
plot(sphere.collect2[[35]][["diag"]])
lines(density(sphe2.per), xlim=c(0, 3), ylim=c(0,100))
sphe2.per.den<-density(sphe2.per)
plot(density(sphe2.per))


detach(package:TDA)
require(ks)
sphe2.per.kde<-kde(sort(sphe2.per))
plot(sphe2.per.kde)
ran<-seq(min(sphe2.per),max(sphe2.per),0.01)
sphe2.per.pl<-plot(ran, pkde(fhat = sphe2.per.kde, q = ran), type="l",xlab="Value",ylab="p")
sphe2.per.thre<-qkde(0.9, sphe2.per.kde)
abline(v=sphe2.per.thre)

sphe2.land<-calcLandscape(list(diagram=sphere.collect2[[35]][["diag"]]))
sphe2.per.thre2<-threshPerKDE((list(diagram=sphere.collect2[[35]][["diag"]])), 0.9)
abline(h=sphe2.per.thre/2, col=6)
abline(h=sphe2.per.thre/4, col=6)

thre.test<-threshPerKDE(list(diagram=sphere.collect2[[35]][["diag"]]))


#KDEによる閾値設定でベッチ数推定試し
#関数がうまくいくか試す用の球状データセット
sphere.collect3<- lapply(1:3, function(i){
  sphere <- sphereUnif(400, 2, 1)
  #cat("nsample=", nsample, "data", i, "\n")
  return(list(nsample = 400, noizyX = sphere, diag = 0))
})
save(sphere.collect3, file = "./data/sphere.collect3")

spre3.aggr<-proposedMethodOnly.kde(sphere.collect3, 2, 3, 10, 300)

data.collect1.land<-calcLandscape(list(diagram=data.collect[[1]][["diag"]]))
dacole.thre<-threshPerKDE((list(diagram=data.collect[[1]][["diag"]])), 0.9)
abline(h=dacole.thre/2, col=6)
abline(h=dacole.thre/4, col=6)

torus10.aggr.test<-proposedMethodOnly.kde(torus.collect10[1], 2, 3, 10)
data1.aggr<-proposedMethodOnly.kde(data.collect, 2, 3, 10)
data1.dim1<-plotAggr(aggr = data1.aggr, collect = data.collect, correct = 2, dim = 1, compare = F, capture = "proposed+kde")
data1.dim2<-plotAggr(aggr = data1.aggr, collect = data.collect, correct = 1, dim = 2, compare = F, capture = "proposed+kde")
data1.holes1<-cyclenumber(data1.aggr[[1]], compare = F)
data1.holes2<-cyclenumber(data1.aggr[[2]], compare = F)

thre.tst<-mapply(function(k)threshPerKDE(data.collect[[k]][["diag"]]), 1:10)
thresh.test2<-threshPerKDE(data.collect[[3]][["diag"]])
data3per<-calcPerdim1dim2(data.collect[[3]][["diag"]])
data3per.kde<-kde(sort(data3per))
ran2<-seq(min(data3per),max(data3per),0.01)
data3per.pl<-plot(ran, pkde(fhat = data3per.kde, q = ran), type="l",xlab="Value",ylab="p")
data3.thre<-qkde(p = 0.9, fhat = data3per.kde)
data3.thre<-qkdeSorted(p = 0.9, fhat = data3per.kde)
data3.land<-calcLandscape(list(diagram=data.collect[[3]][["diag"]]))
abline(h=data3.thre/2, col=6)
abline(h=data3.thre/4, col=6)

sphe2.land26<-calcLandscape(list(diagram=sphere.collect2[[36]][["diag"]]))
sphe2.thre26<-threshPerKDE(list(diagram=sphere.collect2[[36]][["diag"]]), 0.9)
abline(h=sphe2.thre26/2, col=6)
abline(h=sphe2.thre26/4, col=6)
sphe2.nsample<-sapply(sphere.collect2, function(sphere){return(sphere[["nsample"]])})

sphe2.land99<-calcLandscape(list(diagram=sphere.collect2[[99]][["diag"]]))
sphe2.per99<-calcPerdim1dim2(sphere.collect2[[99]][["diag"]])
hist(sphe2.per99, breaks=seq(0,3,0.05), main="Histogram of sphere hole", xlab="persistent range", ylim=c(0,50), col="#993435")
sphe2.per99.dim1<-calcper(sphere.collect2[[99]][["diag"]],1)
hist(sphe2.per99.dim1, breaks=seq(0,3,0.05), main="Histogram of sphere hole", xlab="persistent range", ylim=c(0,50), col="#993435")
per99.mean<-mean(sphe2.per99)
abline(v=per99.mean*2)
abline(v=per99.mean)
sphe2.thre99<-threshPerKDE(list(diagram=sphere.collect2[[99]][["diag"]]), 0.9)
abline(v=sphe2.thre99, col=4)
sphe2.per99.kde<-ks::kde(sort(sphe2.per99))
plot(sphe2.per99.kde)

torus2.land2<-calcLandscape(torus.aggr2[["Xdiag"]][[1]])

#分散を使えないか試し
data3per<-calcPerdim1dim2(data.collect[[3]][["diag"]])
data3per.kde<-ks::kde(sort(log(data3per)))
plot(data3per.kde)
data3per.kde2<-ks::kde(sort(data3per))
plot(data3per.kde2)
data3per.logthre<-threshPerKDE(list(diagram=data.collect[[3]][["diag"]]), 0.9)
abline(v=data3per.logthre)
data3per.thre<-exp(data3per.logthre)
data.collect3.land<-calcLandscape(list(diagram=data.collect[[3]][["diag"]]))
abline(h=data3per.thre/2, col=6)
abline(h=data3per.thre/4, col=6)
data3.per<-calcPerdim1dim2(data.collect[[3]][["diag"]])
hist(data3.per, breaks=seq(0,3,0.05), main="Histogram of noizy torus hole", xlab="persistent range", ylim=c(0,50), col="#993435")

data3.per.mean<-mean(data3.per)
data3.per.sd<-sd(data3.per)
data3.per.scaled<-scale(data3.per)
hist(data3.per.scaled, breaks=seq(-1, 7, 0.1), main="Histogram of noizy torus hole", xlab="persistent range", ylim=c(0,50), col="#993435")

data3.land<-calcLandscape(list(diagram=data.collect[[3]][["diag"]]))
abline(h=data3.thre/2, col=6)
abline(h=data3.thre/4, col=6)
abline(h=data3.thre/2, col=6)
abline(h=(data3.per.mean+(data3.per.sd*2))/4, col=6)
abline(h=(data3.per.mean+(data3.per.sd*2))/2, col=6)

sphe2.per26<-calcPerdim1dim2(sphehole2.aggr[["Xdiag"]][[26]])
sphe2.land26<-calcLandscape(sphehole2.aggr[["Xdiag"]][[26]])
sphe2.per26.mean<-mean(sphe2.per26)
sphe2.per26.sd<-sd(sphe2.per26)
abline(h=(sphe2.per26.mean+sphe2.per26.sd*2)/2, col=6)
abline(h=(sphe2.per26.mean+sphe2.per26.sd*2)/4, col=6)
hist(sphe2.per26, breaks=seq(0,3,0.05), main="Histogram of sphere hole", xlab="persistent range", ylim=c(0,100), col="#993435")

data3.per.quan<-quantile(data3.per)
sphe2.per26.quan<-quantile(sphe2.per26)
plot(sphehole2.aggr[["Xdiag"]][[26]][[1]])
sphe2.per26.dim0<-calcper(sphehole2.aggr[["Xdiag"]][[26]], 0)
sphe2.per26.dim1<-calcper(sphehole2.aggr[["Xdiag"]][[26]], 1)
sphe2.per26.dim2<-calcper(sphehole2.aggr[["Xdiag"]][[26]], 2)

data3.per.dim0<-calcper(data.collect[[3]][["diag"]], 0)
data3.per.dim1<-calcper(data.collect[[3]][["diag"]], 1)
data3.per.dim2<-calcper(data.collect[[3]][["diag"]], 2)
