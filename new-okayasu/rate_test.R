#データ密度と推定制度の関係を調べる
#とりあえず500点トーラスから

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

torus.collect15<- lapply(1:100, function(i){
  nsample <- 500
  #var <- runif(1, var.min, var.max)
  #noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, noizyX = torus, diag = 0))
})
save(torus.collect15, file = "./data/torus.collect15")

torus15.aggr<-proposedMethodOnly(torus.collect15, 2, 3, 10)

save(torus.collect15, file = "./data/torus15.aggr")


torus.collect15.400<-subsampleExclude(torus.collect15, 400)

figurePlot(torus.collect15.400[[1]][["noizyX"]])

#データ密度とベッチ数推定成功率を見る関数づくり
#単位円で試し
circle.test<-circleUnif(500, 1)
plot(circle.test)
circle.test.diag<-ripsDiag(circle.test, 1, 2)
plot(circle.test.diag[[1]])

circle.tst.collect<- lapply(1:3, function(i){
  nsample <- 500
  #var <- runif(1, var.min, var.max)
  #noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
  torus <- circleUnif(500, 1)
  return(list(nsample = nsample, noizyX = torus, diag = 0))
})

circle.tst.aggrlist<-successTimes(circle.tst.collect, 200, 2, 1, 2)

circle.tst.aggr<-proposedMethodOnly.mk1(circle.tst.collect, 1, 2, 10)

circle.tst.subcollect.300<-subsampleExclude(circle.tst.collect, 300)
circle.tst.subaggr<-proposedMethodOnly.mk1(circle.tst.subcollect.300, 1, 2, 10)

circle.tst.suc300<-successRates(circle.tst.aggrlist, 1, 1)

per.set<-lapply(data.collect[1:10], function(data){
  return(calcPerdim1dim2(data[["diag"]]))
})

circle.tst.rates<-succesCheck(circle.tst.collect, nsubset=200, times=2, maxdim=1, mxscale=2, corrects=1)



test.sub<-subsampleExclude(torus15.300sub, nsub = 100)
save(test.sub, file="./data/test_subA.RData")
figurePlot(test.sub[[1]][["noizyX"]])
save2File(test.sub)


torus15.dim1<-cyclenumber(torus15.aggr[[1]], compare=F)
torus15.dim2<-cyclenumber(torus15.aggr[[2]], compare=F)

torus15.sucs<-successRates(list(append(torus15.aggr, list(dim1cycles=torus15.dim1, dim2cycles=torus15.dim2))), 1, 2)
torus15.sucsdim2<-successRates(list(append(torus15.aggr, list(dim1cycles=torus15.dim1, dim2cycles=torus15.dim2))), 2, 1)


torus15.rates<-succesCheck(torus.collect15, nsubset=c(200, 250, 300, 350, 400, 450),
                    times=5, maxdim=2, mxscale=3, corrects=c(2,1), origin = T)

#データ点数と成功率の関係を見る関数の不備を直す
trus200.tst<- lapply(1:2, function(i){
  nsample <- 200
  #var <- runif(1, var.min, var.max)
  #noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, noizyX = torus, diag = 0))
})
figurePlot(trus200.tst[[1]][["noizyX"]])

trus200.rates<-succesCheck(trus200.tst, nsubset=100,
                           times=2, maxdim=2, mxscale=3, corrects=c(2,1), origin = F)

trus200.diag<-ripsDiag(trus200.tst[[1]][["noizyX"]], 2, 3)

trus200.list<-list(trus200.tst[[1]][["noizyX"]])
class(trus200.list)<-"bootsSamples"
trus200.boot<-bootstrap.homology.mk3(trus200.B, 2, 3)

trus200.B<-bootstrapper(trus200.tst[[1]][["noizyX"]],100,10)
trus200.boot<-bootstrap.homology.mk3(trus200.B, 2, 3)

trus200.time<-successTimes(trus200.tst, 100, 2, 2, 3)
trus200.aggr<-proposedMethodOnly.mk1(trus200.tst, 2, 3, 3)
trus200.rate<-successRates(trus200.time, 1, 2)
trus200.rate2<-successRates(trus200.time, 2, 1)

trus200.intime<-interSuccessTimes(trus200.tst, 100, 1, 2, 3)

trus200.inrate<-successRates(trus200.intime, dim = 1, correct = 2, inter=T)

trus200.inrates<-interpoCheck(trus200.tst, nsubset=100,
                              times=2, maxdim=2, mxscale=3, corrects=c(2,1), origin = F)

trus200tst.subs<-lapply(1:2, function(k)subsampleExclude(trus200.tst, 100))

trus200tst.aggrs<-lapply(1:2, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(trus200tst.subs[[k]], 2, 3, 10))})

trus200tst.cycs<-lapply(trus200tst.aggrs, function(aggr){
  dim1<-cyclenumber(aggr[[1]], compare = F)
  dim2<-cyclenumber(aggr[[2]], compare = F)
  return(list(dim1=dim1, dim2=dim2))})

trus200tst.rates<-aggrSuccessRates(trus200tst.aggrs, correct = c(2,1))
#######################################
torus15.300rates<-succesCheck(torus.collect15, nsubset=300,
                           times=5, maxdim=2, mxscale=3, corrects=c(2,1))

trus15.diag1<-ripsDiag(torus.collect15[[1]][["noizyX"]], 2, 3)
par(cex.lab=1.6, cex.axis=1.6)
plot(trus15.diag1[[1]])
save2Rdata(trus15.diag1)

trus15.diag1a<-compute_pd(torus.collect15[[1]][["noizyX"]], 2, 3)
autoplot(trus15.diag1a)

trus15.pl1<-compute_pl(trus15.diag1a)
autoplot(trus15.pl1)

trus15.smpl1<-compute_smooth_pl(trus15.pl1)
autoplot(trus15.smpl1)

trus15.betti1<-count_smooth_maximal(trus15.pl1)

trus15.pL1<-calcLandscape(trus15.diag1)
plot(trus15.pL1$tseq, trus15.pL1$Land.dim1, type = "l", col=2, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 1), cex.lab=1.6, cex.axis=1.6, lwd=2)
par(new=T, mgp=c(2.5, 1, 0))
plot(trus15.pL1$tseq, trus15.pL1$Land.dim2, type = "l", col=3, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 1), cex.lab=1.6, cex.axis=1.6, lwd=2)
par(mgp=c(2.5, 1, 0))
legend("topleft", legend = c("dim1", "dim2"), col = c(2, 3), lty = 1, cex=1.5, lwd = 2)

peak<-calc.landscape.peak(trus15.pL1$Land.dim1, thresh = calcPerMeanNoDouble(trus15.diag1[[1]]), tseq = trus15.pL1$tseq, show = T)

figurePlot(torus.collect15[[1]][["noizyX"]])
rgl.postscript("./data/torus.eps", fmt="eps" ) 

torus15.350rates<-succesCheck(torus.collect15, nsubset=350,
                              times=5, maxdim=2, mxscale=3, corrects=c(2,1))

torus15.340rates<-succesCheck(torus.collect15, nsubset=340,
                              times=5, maxdim=2, mxscale=3, corrects=c(2,1))

torus15.330rates<-succesCheck(torus.collect15, nsubset=330,
                              times=5, maxdim=2, mxscale=3, corrects=c(2,1))

torus15.320rates<-succesCheck(torus.collect15, nsubset=320,
                              times=5, maxdim=2, mxscale=3, corrects=c(2,1))

save(torus15.300rates, file="./data/trus15_300rates.RData")
save(torus15.350rates, file="./data/trus15_350rates.RData")
save(torus15.340rates, file="./data/trus15_340rates.RData")
save(torus15.330rates, file="./data/trus15_330rates.RData")
save(torus15.320rates, file="./data/trus15_320rates.RData")

sucrate300sub.dim2<-sapply(torus15.300subs.rate, function(rate){return(rate[2])})
sucrate320sub.dim2<-sapply(torus15.320rates[["320sub"]][["dim2"]], function(rate){return(rate[1])})
sucrate350sub.dim2<-sapply(torus15.350rates[["350sub"]][["dim2"]], function(rate){return(rate[1])})
sucrate330sub.dim2<-sapply(torus15.330rates[["330sub"]][["dim2"]], function(rate){return(rate)})
sucrate340sub.dim2<-sapply(torus15.340rates[["340sub"]][["dim2"]], function(rate){return(rate)})

sucrate310sub.dim2<-sapply(torus15.310subs.rates, function(rate){return(rate[["dim2rate"]])})
sucrate300sub2.dim2<-sapply(torus15.300subs.rate, function(rate){return(rate[["dim2rate"]])})

#精度のプロット
suctrate.dim2<-list("300"=sucrate300sub2.dim2,
                    "310"=sucrate310sub.dim2,
                           "320"=sucrate320sub.dim2,
                           "330"=sucrate330sub.dim2,
                           "340"=sucrate340sub.dim2,
                           "350"=sucrate350sub.dim2)

suctrate.dim2[["300"]]<-unlist(sucrate300sub.dim2)
suctrate.dim2[["310"]]<-unlist(sucrate310sub.dim2)

suctrate.dim2 %>% bind_cols() %>% gather(data, value) %>% ggplot(aes(data, value)) + geom_violin() + geom_point()

save2Rdata(suctrate.dim2)

sucrate.dim2.tidy<-suctrate.dim2 %>% bind_cols() %>% gather(data, value)

rate.dim2<-sucrate.dim2.tidy[,"value"]

plot(sucrate.dim2.tidy, pch=16, cex.axis=1.6, xlab="Data Points", ylab="Success Rates", cex.lab=1.6)

sucdim2.mean<-sapply(suctrate.dim2, function(rate)mean(rate))
lines(seq(300, 350, by=10), sucdim2.mean)

sucdim2.sd<-sapply(suctrate.dim2, function(rate)sd(rate))
lines(seq(300, 350, by=10), sucdim2.mean-sucdim2.sd, lty="dashed")
lines(seq(300, 350, by=10), sucdim2.mean+sucdim2.sd, lty="dashed")
####################


torus15.300inrates<-interpoCheck(torus.collect15[1:50], nsubset=300,
                              times=1, maxdim=2, mxscale=3, corrects=c(2,1), origin = F)
save(torus15.300inrates, file="./data/trus15_300inrates.RData")

##ボロノイ領域の中心に点を打つ
torus15.300inrates.cen<-interpoCheck(torus.collect15[1:50], nsubset=300,
                                 times=1, maxdim=2, mxscale=3, corrects=c(2,1), origin = F, border = F)

#サブサンプルセットを作る################
torus15.300sub<-subsampleExclude(torus.collect15, nsub = 300)
save(torus15.300sub, file="./data/torus15_300sub.RData")
torus15.310sub<-subsampleExclude(torus.collect15, nsub = 310)
save(torus15.310sub, file="./data/torus15_310sub.RData")
torus15.320sub<-subsampleExclude(torus.collect15, nsub = 320)
save(torus15.320sub, file="./data/torus15_320sub.RData")
torus15.330sub<-subsampleExclude(torus.collect15, nsub = 330)
save(torus15.330sub, file="./data/torus15_330sub.RData")
torus15.340sub<-subsampleExclude(torus.collect15, nsub = 340)
save(torus15.340sub, file="./data/torus15_340sub.RData")
torus15.350sub<-subsampleExclude(torus.collect15, nsub = 350)
save(torus15.350sub, file="./data/torus15_350sub.RData")
torus15.400sub<-subsampleExclude(torus.collect15, nsub = 400)
save(torus15.400sub, file="./data/torus15_400sub.RData")
torus15.450sub<-subsampleExclude(torus.collect15, nsub = 450)
save(torus15.450sub, file="./data/torus15_450sub.RData")
torus15.250sub<-subsampleExclude(torus.collect15, nsub = 250)
save(torus15.250sub, file="./data/torus15_250sub.RData")
torus15.200sub<-subsampleExclude(torus.collect15, nsub = 200)
save(torus15.200sub, file="./data/torus15_200sub.RData")
#############################

torus15.300insub<-torus15.300sub
for (l in 1:length(torus15.300sub)) {
  inter.oricord<-voronoiInterpo(torus15.300sub[[l]][["noizyX"]], 15, border = F)
  torus15.300insub[[l]][["noizyX"]]<-rbind(torus15.300insub[[l]][["noizyX"]], inter.oricord)
  torus15.300insub[[l]][["nsample"]]<-nrow(torus15.300insub[[l]][["noizyX"]])
  debugText(l, torus15.300insub[[l]][["nsample"]])
}

#サブサンプルリストを作る2############
torus15.310subs<-lapply(1:5, function(k)subsampleExclude(torus.collect15, nsub = 310))
save(torus15.310subs, file="./data/torus15_310subs.RData")

torus15.300subs<-lapply(1:5, function(k)subsampleExclude(torus.collect15, nsub = 300))
save(torus15.300subs, file="./data/torus15_300subs.RData")

######################################

torus15.310subs.aggrs<-lapply(1:length(torus15.310subs), function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly.mk1(torus15.310subs[[k]], 2, 3, 10))})

save(torus15.310subs.aggrs, file="./data/torus15_310subs_aggrs.RData")

torus15.310subs.rates<-aggrSuccessRates(torus15.310subs.aggrs, correct = c(2,1))
save(torus15.310subs.rates, file="./data/torus15_310subs_rates.RData")

torus15.300subs.aggrs<-lapply(1:length(torus15.300subs), function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly.mk1(torus15.300subs[[k]], 2, 3, 10))})

save(torus15.300subs.aggrs, file="./data/torus15_300subs_aggrs.RData")


torus15.300subs.rates<-aggrSuccessRates(torus15.300subs.aggrs, correct = c(2,1))
save(torus15.300subs.rates, file="./data/torus15_300subs_rates.RData")

#サブサンプルに補間
torus15.300insub1<-intering(torus15.300subs[[1]])
save(torus15.300insub1, file="./data/torus15_300insub1.RData")
##############################

#補間したサブサンプルを分析
torus15.300insub1.aggr<-proposedMethodOnly(torus15.300insub1, 2, 3, 10)
save(torus15.300insub1.aggr, file="./data/torus15_300insub1_aggr.RData")

torus15.300insubs1.rate<-aggrSuccessRates(list(torus15.300insub1.aggr), correct = c(2,1))


torus15.300insubs<-lapply(torus15.300subs, function(sub)intering(sub))
save(torus15.300insubs, file="./data/torus15_300insubs.RData")

torus15.300insubs2_3.aggrs<-lapply(2:3, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.300insubs[[k]], 2, 3, 10))})

save(torus15.300insubs2_3.aggrs, file="./data/torus15_300insubs2_3_aggrs.RData")

torus15.300insubs2_3.rate<-aggrSuccessRates(torus15.300subs2_3.aggrs, correct = c(2,1))



