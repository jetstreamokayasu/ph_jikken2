#データ密度と推定制度の関係を調べる
#とりあえず500点トーラスから
#ごちゃごちゃになってしまったので整理

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
library(seephacm)


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

#サブサンプルリストを作る3############
torus15.310subs<-lapply(1:5, function(k)subsampleExclude(torus.collect15, nsub = 310))
save(torus15.310subs, file="./data/torus15_310subs.RData")

torus15.300subs<-lapply(1:5, function(k)subsampleExclude(torus.collect15, nsub = 300))
save(torus15.300subs, file="./data/torus15_300subs.RData")

######################################

#サブサンプルを分析##############
torus15.300subs.aggrs<-lapply(1:length(torus15.300subs), function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.300subs[[k]], 2, 3, 10))})
save(torus15.300subs.aggrs, file="./data/torus15_300subs_aggrs.RData")

torus15.300subs.rate<-aggrSuccessRates(torus15.300subs.aggrs, correct = c(2,1))

###############################

#サブサンプルに補間
torus15.300insubs<-lapply(torus15.300subs, function(sub)intering(sub))
save(torus15.300insubs, file="./data/torus15_300insubs.RData")

torus15.310insubs<-lapply(torus15.310subs, function(sub)intering(sub))
save2Rdata(torus15.310insubs)
##############################

#補間したサブサンプルを分析
torus15.300subs4_5.aggrs<-lapply(4:5, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.300insubs[[k]], 2, 3, 10))})

save(torus15.300subs4_5.aggrs, file="./data/torus15_300subs4_5_aggrs.RData")

torus15.300insubs4_5rate<-aggrSuccessRates(torus15.300subs4_5.aggrs, correct = c(2,1))

#凸包外にある補間点を除く補間
torus15.300insubs1_3<-lapply(torus15.300subs[1:3], function(sub)intering(sub))
save2Rdata(torus15.300insubs1_3)

torus15.300sub1_27<-voronoiInterpo(torus15.300subs[[1]][[27]][["noizyX"]], 15)
torus15.300sub2_56<-voronoiInterpo(torus15.300subs[[2]][[56]][["noizyX"]], 15)

torus15.300insubs1_3.aggrs<-lapply(1:3, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.300insubs1_3[[k]], 2, 3, 10))})

torus15.300sub1_36<-proposedMethodOnly(torus15.300insubs1_3[[1]][36], maxdim = 2, maxscale = 3, samples = 10)

torus15.300insubs2.aggr<-proposedMethodOnly(torus15.300insubs1_3[[2]], maxdim = 2, maxscale = 3, samples = 10)
save2Rdata(torus15.300insubs2.aggr)
torus15.300insubs2.rate<-aggrSuccessRates(list(torus15.300insubs2.aggr), correct=c(2,1))

torus15.300insubs1.aggr<-proposedMethodOnly(torus15.300insubs1_3[[1]], maxdim = 2, maxscale = 3, samples = 10)
save2Rdata(torus15.300insubs1.aggr)
torus15.300insubs1.rate<-aggrSuccessRates(list(torus15.300insubs1.aggr), c(2,1))

torus15.300insubs3.aggr<-proposedMethodOnly(torus15.300insubs1_3[[3]], maxdim = 2, maxscale = 3, samples = 10)
save2Rdata(torus15.300insubs3.aggr)
torus15.300insubs3.rate<-aggrSuccessRates(list(torus15.300insubs3.aggr), c(2,1))

#300点
torus15.300insubs<-lapply(torus15.300subs, function(sub)intering(sub))
save2Rdata(torus15.300insubs)

#310点
torus15.310insubs1_3.aggrs<-lapply(1:3, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.310insubs[[k]], 2, 3, 10))})
save2Rdata(torus15.310insubs1_3.aggrs)

torus15.310insubs1_3.rates<-aggrSuccessRates(torus15.310insubs1_3.aggrs, c(2,1))

torus15.310insubs4_5.aggrs<-lapply(4:5, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.310insubs[[k]], 2, 3, 10))})
save2Rdata(torus15.310insubs4_5.aggrs)

torus15.310insubs4_5.rates<-aggrSuccessRates(torus15.310insubs4_5.aggrs, c(2,1))

#320点
torus15.320subs<-lapply(1:5, function(k)subsampleExclude(torus.collect15, nsub = 320))
save(torus15.320subs, file="./data/torus15_320subs.RData")

torus15.320insubs<-lapply(torus15.320subs, function(sub)intering(sub))
save2Rdata(torus15.320insubs)

torus15.320insubs1_3.aggrs<-lapply(1:3, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.320insubs[[k]], 2, 3, 10))})
save2Rdata(torus15.320insubs1_3.aggrs)

torus15.320insubs1_3.rates<-aggrSuccessRates(torus15.320insubs1_3.aggrs, c(2,1))

torus15.320insubs4_5.aggrs<-lapply(4:5, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.320insubs[[k]], 2, 3, 10))})
save2Rdata(torus15.320insubs4_5.aggrs)

torus15.320insubs4_5.rates<-aggrSuccessRates(torus15.320insubs4_5.aggrs, c(2,1))

#330点
torus15.330subs<-lapply(1:5, function(k)subsampleExclude(torus.collect15, nsub = 330))
save2Rdata(torus15.330subs)

torus15.330insubs<-lapply(torus15.330subs, function(sub)intering(sub))
save2Rdata(torus15.330insubs)

torus15.330insubs1_3.aggrs<-lapply(1:3, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.330insubs[[k]], 2, 3, 10))})
save2Rdata(torus15.330insubs1_3.aggrs)

torus15.330insubs1_3.rates<-aggrSuccessRates(torus15.330insubs1_3.aggrs, c(2,1))

torus15.330insubs4_5.aggrs<-lapply(4:5, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.330insubs[[k]], 2, 3, 10))})
save2Rdata(torus15.330insubs4_5.aggrs)
save(torus15.330insubs4_5.aggrs, file="torus15.330insubs4_5.aggrs.RData")

torus15.330insubs4_5.rates<-aggrSuccessRates(torus15.330insubs4_5.aggrs, c(2,1))

#340点
torus15.340subs<-lapply(1:5, function(k)subsampleExclude(torus.collect15, nsub = 340))
save2Rdata(torus15.340subs)

torus15.340insubs<-lapply(torus15.340subs, function(sub)intering(sub))
save2Rdata(torus15.340insubs)

torus15.340insubs1_3.aggrs<-lapply(1:3, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.340insubs[[k]], 2, 3, 10))})
save2Rdata(torus15.340insubs1_3.aggrs)

torus15.340insubs1_3.rates<-aggrSuccessRates(torus15.340insubs1_3.aggrs, c(2,1))

torus15.340insubs4_5.aggrs<-lapply(4:5, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.340insubs[[k]], 2, 3, 10))})
save2Rdata(torus15.340insubs4_5.aggrs)

torus15.340insubs4_5.rates<-aggrSuccessRates(torus15.340insubs4_5.aggrs, c(2,1))

#350点

torus15.350subs<-lapply(1:5, function(k)subsampleExclude(torus.collect15, nsub = 350))
save2Rdata(torus15.350subs)

torus15.350insubs1_3.rates<-aggrSuccessRates(torus15.350insubs1_3.aggrs, c(2,1))

torus15.350insubs4_5.aggrs<-lapply(4:5, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.350insubs[[k]], 2, 3, 10))})
save2Rdata(torus15.350insubs4_5.aggrs)

torus15.350insubs4_5.rates<-aggrSuccessRates(torus15.350insubs4_5.aggrs, c(2,1))

#250点

torus15.250subs<-lapply(1:5, function(k)subsampleExclude(torus.collect15, nsub = 250))
save2Rdata(torus15.250subs)
save(torus15.250subs, file="./data/torus15.250sub.RData")

torus15.250subs1_3.aggrs<-lapply(1:3, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.250subs[[k]], 2, 3, 10))})
save2Rdata(torus15.250subs1_3.aggrs)
save(torus15.250subs1_3.aggrs, file="./data/torus15.250sub1_3.aggrs.RData")

torus15.250subs1_3.rates<-aggrSuccessRates(torus15.250subs1_3.aggrs, c(2,1))

torus15.250insubs<-lapply(torus15.250subs, function(sub)intering(sub))
save2Rdata(torus15.250insubs)

torus15.250insubs1_3.aggrs<-lapply(1:3, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.250insubs[[k]], 2, 3, 10))})
save2Rdata(torus15.250insubs1_3.aggrs)

torus15.250insubs1_3.rates<-aggrSuccessRates(torus15.250insubs1_3.aggrs, c(2,1))

##精度のプロット2
suctrate.dim2_2<-list("300"=unlist(sucrate300sub.dim2),
                      "310"=sucrate310sub.dim2,
                      "320"=suctrate.dim2[["320"]],
                      "330"=suctrate.dim2[["330"]],
                      "340"=suctrate.dim2[["340"]],
                      "350"=suctrate.dim2[["350"]])


torus15.350insubs<-lapply(torus15.350subs, function(sub)intering(sub))
save2Rdata(torus15.350insubs)

torus15.350insubs1_3.aggrs<-lapply(1:3, function(k){
  cat("list", k, "calc\n")
  return(proposedMethodOnly(torus15.350insubs[[k]], 2, 3, 10))})
save2Rdata(torus15.350insubs1_3.aggrs)


torus15.350insubs1_3.rates<-aggrSuccessRates(torus15.350insubs1_3.aggrs, c(2,1))

suctrate.dim2_2 %>% bind_cols() %>% gather(data, value) %>% ggplot(aes(data, value)) + geom_violin() + geom_point()

sucrate.dim2.tidy_2<-suctrate.dim2_2 %>% bind_cols() %>% gather(data, value)

rate.dim2<-sucrate.dim2.tidy_2[,"value"]

oldpar <- par(no.readonly=T)

plot(sucrate.dim2.tidy_2, pch=16, cex.axis=1.6, xlab="Data Density", ylab="Success Rates", cex.lab=1.6, ylim=c(0.2, 1.0), xaxt="n")
par(family = "serif")
par(family="")
axis(side=1, at=seq(300, 350, by=10), labels=c(paste0(seq(30, 35), "/(pi^2)")), cex.axis=1.1)

sucdim2.mean_2<-sapply(suctrate.dim2_2, function(rate)mean(rate))
lines(seq(300, 350, by=10), sucdim2.mean_2)

sucdim2.sd_2<-sapply(suctrate.dim2_2, function(rate)sd(rate))
lines(seq(300, 350, by=10), sucdim2.mean_2-sucdim2.sd_2, lty="dashed")
lines(seq(300, 350, by=10), sucdim2.mean_2+sucdim2.sd_2, lty="dashed")

#補間後の成功率まとめ
insub300.rate<-c(torus15.300insubs1.rate[[1]][["dim2rate"]], torus15.300insubs2.rate[[1]][["dim2rate"]], torus15.300insubs3.rate[[1]][["dim2rate"]], 
                 torus15.300insubs4_5rate[[1]][["dim2rate"]], torus15.300insubs4_5rate[[2]][["dim2rate"]])

insub310.rate<-sapply(append(torus15.310insubs1_3.rates, torus15.310insubs4_5.rates), function(rate){return(rate[[2]])})

insub320.rate<-sapply(append(torus15.320insubs1_3.rates, torus15.320insubs4_5.rates), function(rate){return(rate[[2]])})

insub330.rate<-sapply(append(torus15.330insubs1_3.rates, torus15.330insubs4_5.rates), function(rate){return(rate[[2]])})

insub340.rate<-sapply(append(torus15.340insubs1_3.rates, torus15.340insubs4_5.rates), function(rate){return(rate[[2]])})

insub350.rate<-sapply(append(torus15.350insubs1_3.rates, torus15.350insubs4_5.rates), function(rate){return(rate[[2]])})

insub.rates<-list("300"=insub300.rate,
                  "310"=insub310.rate,
                  "320"=insub320.rate,
                  "330"=insub330.rate,
                  "340"=insub340.rate,
                  "350"=insub350.rate)

save2Rdata(insub.rates)

points(rep(300, 5), insub300.rate, col=2, pch=16)
points(rep(310, 5), insub310.rate, col=2, pch=16)
points(rep(320, 5), insub320.rate, col=2, pch=16)
points(rep(330, 5), insub330.rate, col=2, pch=16)
points(rep(340, 5), insub340.rate, col=2, pch=16)
points(rep(350, 5), insub350.rate, col=2, pch=16)

insubdim2.mean<-sapply(insub.rates, function(rate)mean(rate))
lines(seq(300, 350, by=10), insubdim2.mean, col=2)

insubdim2.sd<-sapply(insub.rates, function(rate)sd(rate))
lines(seq(300, 350, by=10), insubdim2.mean-insubdim2.sd, lty="dashed", col=2)
lines(seq(300, 350, by=10), insubdim2.mean+insubdim2.sd, lty="dashed", col=2)

#補間プログラムのパッケージ実験
torus.collect16<- lapply(1:100, function(i){
  nsample <- 500
  #var <- runif(1, var.min, var.max)
  #noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, X = torus, diag = 0))
})
save2Rdata(torus.collect16)

trs16_vic1<-get_vicinity(trs16_dist, 1, 15)
trs16_dist<-dist(torus.collect16[[1]][[2]])
figurePlot3d(torus.collect16[[1]][[2]][-trs16_vic1, ])
points3d(torus.collect16[[1]][[2]][trs16_vic1, ], col=3)
