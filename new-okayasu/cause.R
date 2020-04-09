require(phacm)
require(tidyverse)
require(TDA)
require(myfs)
require(rgl)
library(seephacm)
library(tictoc)
library(devtools)
library(interpo3d)

#なぜ1次ベッチ数の精度が下がるのか調べる
#350点トーラスで比較
#torus.collect18とtorus_colle_set1が350点の補間前トーラス

#補間前に正解なのと補間後に正解しなかったデータセット比較
crrct1_10<-which(torus_colset1_aggrs[[1]][[1]] >= 1.5 & torus_colset1_aggrs[[1]][[1]] < 2.5)

inwrng1_10<-which(torus350_incolle_aggrs[[2]][[1]] < 1.5 | torus_colset1_aggrs[[1]][[1]] >= 2.5)
#上のwhichの中身間違えている可能性大

c_inw_1_10<-intersect(crrct1_10, inwrng1_10)[1:10]

trs350_2_cw1_10_pd<-lapply(c_inw_1_10, function(k)ripsDiag(torus350_colle_set[[2]][[k]][["noizyX"]], 2, 3, printProgress = T))
plotPDs(trs350_2_cw1_10_pd)

trs350_2_cw1_10_pls<-lapply(1:10, function(k)calcLandscape(trs350_2_cw1_10_pd[[k]]))
plot_lands(trs350_2_cw1_10_pls, 1)

trs_in350_2_cw1_10_pd<-lapply(c_inw_1_10, function(k)ripsDiag(torus350_incolle_set[[2]][[k]][["noizyX"]], 2, 3, printProgress = T))
plotPDs(trs_in350_2_cw1_10_pd)

trs_in350_2_cw1_10_pls<-lapply(1:10, function(k)calcLandscape(trs_in350_2_cw1_10_pd[[k]]))
plot_lands(trs_in350_2_cw1_10_pls, 1)


#新手法の補間の誤差を調べる
trs_in350_2_26_der<-torus_disterror(torus350_incolle_set[[2]][[26]][["noizyX"]], maxr = 2.5, minr = 1, nps = 350)
hist(trs_in350_2_26_der, col="#993435")

oldpar <- par(no.readonly = TRUE)  

trs_in350_1_w1_10_ders<-lapply(inwrng1_10, function(i)torus_disterror(torus350_incolle_set[[2]][[i]][["noizyX"]], maxr = 2.5, minr = 1, nps = 350))
par(mgp=c(2.5,1,0))
boxplot(trs_in350_1_w1_10_ders, xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6)


trs_in300_1_1_10_ders<-lapply(1:10, function(i)torus_disterror(torus300_1_1[[i]][["noizyX"]], maxr = 2.5, minr = 1, nps = 300))
par(mgp=c(2.5,1,0))
boxplot(trs_in300_1_1_10_ders, xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6)


#サブサンプルのPDを確かめる
##補間前
trs_350subs2_26_subs<-lapply(1:10, function(k){
  data<-torus_colle_set1[[1]][[26]][["noizyX"]][sample(350, 350*0.8),]
  return(data)
})

trs_350subs2_26_subs_pds<-lapply(1:10, function(k)
  ripsDiag(trs_350subs2_26_subs[[k]], maxdimension = 2, maxscale = 3, printProgress = T))

plotPDs(trs_350subs2_26_subs_pds)

trs_350subs2_26_subs_pls<-lapply(1:10, function(k)calcLandscape(trs_350subs2_26_subs_pds[[k]]))
plot_lands(trs_350subs2_26_subs_pls, dim = 1)


##補間後
trs_in350subs2_26_subs<-lapply(1:10, function(k){
  npoints<-nrow(torus350_incolle_set[[2]][[26]][["noizyX"]])
  data<-torus350_incolle_set[[2]][[26]][["noizyX"]][sample(npoints, npoints*0.8),]
  return(data)
})

trs_in350subs2_26_subs_pds<-lapply(1:10, function(k)
  ripsDiag(trs_in350subs2_26_subs[[k]], maxdimension = 2, maxscale = 3, printProgress = T))

plotPDs(trs_in350subs2_26_subs_pds)

trs_in350subs2_26_subs_pls<-lapply(1:10, function(k)calcLandscape(trs_in350subs2_26_subs_pds[[k]]))
plot_lands(trs_in350subs2_26_subs_pls, dim = 1)


#補間点数の差による計算時間差を調べる
#1点、3点、6点は終了
#新手法で。補間点は遠い順
##300点トーラスで試す
t300_a1_intime<-system.time(torus300_incolle_1_a1<-variable_interpo(torus300_colle_set[[1]], 15, 1))

t300_a2_intime<-system.time(torus300_incolle_1_a2<-variable_interpo(torus300_colle_set[[1]], 15, 2))

t300_a3_intime<-system.time(torus300_incolle_1_a3<-variable_interpo(torus300_colle_set[[1]], 15, 3))

t300_a4_intime<-system.time(torus300_incolle_1_a4<-variable_interpo(torus300_colle_set[[1]], 15, 4))

t300_a5_intime<-system.time(torus300_incolle_1_a5<-variable_interpo(torus300_colle_set[[1]], 15, 5))

#1点、3点、6点は推定終了
t300_a6_intime<-system.time(torus300_incolle_1_a6<-variable_interpo(torus300_colle_set[[1]], 15, 6))

t300_a7_intime<-system.time(torus300_incolle_1_a7<-variable_interpo(torus300_colle_set[[1]], 15, 7))

t300_a8_intime<-system.time(torus300_incolle_1_a8<-variable_interpo(torus300_colle_set[[1]], 15, 8))

torus300_1_test<-all_interpolate(torus300_colle_set[[1]], nvic = 15)

torus300_incolle_1_a6_2<-variable_interpo(torus300_colle_set[[1]], 15, 6)

torus300_incolle_1_a6_2_77<-variable_interpo(list(torus300_colle_set[[1]][[77]]), 15, 6)

intt300_1_a3_time<-system.time(trs300_incolle1_a3_aggr<-proposedMethodOnly(torus300_incolle_1_a3, 2, 3, 10))

intt300_1_a6_time<-system.time(trs300_incolle1_a6_aggr<-proposedMethodOnly(torus300_incolle_1_a6, 2, 3, 10))

intt300_1_a1_time<-system.time(trs300_incolle1_a1_aggr<-proposedMethodOnly(torus300_incolle_1_a1, 2, 3, 10))

#2,4,8点を推定
intt300_1_a8_time<-system.time(trs300_incolle1_a8_aggr<-proposedMethodOnly(torus300_incolle_1_a8, 2, 3, 10))

intt300_1_a2_time<-system.time(trs300_incolle1_a2_aggr<-proposedMethodOnly(torus300_incolle_1_a2, 2, 3, 10))

intt300_1_a4_time<-system.time(trs300_incolle1_a4_aggr<-proposedMethodOnly(torus300_incolle_1_a4, 2, 3, 10))

#1,2,4,6点を2セット推定
##1点追加
torus300_incolle_a1<-lapply(torus300_colle_set, function(k){
  
  time<-system.time(incole<-variable_interpo(k, 15, 1))
  return(list(incole, time=time))
  
})

intt300_1_a1_time2<-system.time(trs300_incolle1_a1_aggr2<-proposedMethodOnly(torus300_incolle_a1[[1]][[1]], 2, 3, 10))

trs300_incole13_a1_aggrs<-lapply(2:3, function(k){
  
  cat("list", k, "calc\n")
  time<-system.time(aggr<-proposedMethodOnly(torus300_incolle_a1[[k]][[1]], 2, 3, 10))
  return(append(aggr, list(time=time)))
  
})
save2Rdata(trs300_incole13_a1_aggrs)


###1点追加4,5セット目推定
{
  trs300_incole45_a1_aggrs<-lapply(4:5, function(k){
    
    cat("list", k, "calc\n")
    time<-system.time(aggr<-proposedMethodOnly(torus300_incolle_a1[[k]][[1]], 2, 3, 10))
    return(append(aggr, list(time=time)))
    
  })
  save2Rdata(trs300_incole45_a1_aggrs)
}

##2点追加
torus300_incolle_a2<-lapply(torus300_colle_set, function(k){
  
  time<-system.time(incole<-variable_interpo(k, 15, 2))
  return(list(incole, time=time))
  
})

###2点追加2,3セット目推定
trs300_incole13_a2_aggrs<-lapply(2:3, function(k){
  
  cat("list", k, "calc\n")
  time<-system.time(aggr<-proposedMethodOnly(torus300_incolle_a2[[k]][[1]], 2, 3, 10))
  return(append(aggr, list(time=time)))
  
})
save2Rdata(trs300_incole13_a2_aggrs)


###2点追加4,5セット目推定
{
trs300_incole45_a2_aggrs<-lapply(4:5, function(k){
  
  cat("list", k, "calc\n")
  time<-system.time(aggr<-proposedMethodOnly(torus300_incolle_a2[[k]][[1]], 2, 3, 10))
  return(append(aggr, list(time=time)))
  
})
save2Rdata(trs300_incole45_a2_aggrs)
}


##3点追加
torus300_incolle_a3<-lapply(torus300_colle_set, function(k){
  
  time<-system.time(incole<-variable_interpo(k, 15, 3))
  return(list(incole, time=time))
  
})

intt300_1_a3_time2<-system.time(trs300_incolle1_a3_aggr2<-proposedMethodOnly(torus300_incolle_a3[[1]][[1]], 2, 3, 10))

###3点追加2,3セット目推定
{
trs300_incole23_a3_aggrs<-lapply(2:3, function(k){
  
  cat("list", k, "calc\n")
  time<-system.time(aggr<-proposedMethodOnly(torus300_incolle_a3[[k]][[1]], 2, 3, 10))
  return(append(aggr, list(time=time)))
  
})
save2Rdata(trs300_incole23_a3_aggrs)
}

###3点追加4,5セット目推定

{
  trs300_incole45_a3_aggrs<-lapply(4:5, function(k){
    
    cat("list", k, "calc\n")
    time<-system.time(aggr<-proposedMethodOnly(torus300_incolle_a3[[k]][[1]], 2, 3, 10))
    return(append(aggr, list(time=time)))
    
  })
  save2Rdata(trs300_incole45_a3_aggrs)
}

##4点追加
torus300_incolle_a4<-lapply(torus300_colle_set, function(k){
  
  time<-system.time(incole<-variable_interpo(k, 15, 4))
  return(list(incole, time=time))
  
})

###4点追加2,3セット目推定
{
  trs300_incole23_a4_aggrs<-lapply(2:3, function(k){
    
    cat("list", k, "calc\n")
    time<-system.time(aggr<-proposedMethodOnly(torus300_incolle_a4[[k]][[1]], 2, 3, 10))
    return(append(aggr, list(time=time)))
    
  })
  save2Rdata(trs300_incole23_a4_aggrs)
}

###4点追加4,5セット目推定

{
  trs300_incole45_a4_aggrs<-lapply(4:5, function(k){
    
    cat("list", k, "calc\n")
    time<-system.time(aggr<-proposedMethodOnly(torus300_incolle_a4[[k]][[1]], 2, 3, 10))
    return(append(aggr, list(time=time)))
    
  })
  save2Rdata(trs300_incole45_a4_aggrs)
}



##6点追加
torus300_incolle_a6<-lapply(torus300_colle_set, function(k){
  
  time<-system.time(incole<-variable_interpo(k, 15, 6))
  return(list(incole, time=time))
  
})


###6点追加2,3セット目推定

{
  trs300_incole23_a6_aggrs<-lapply(2:3, function(k){
    
    cat("list", k, "calc\n")
    time<-system.time(aggr<-proposedMethodOnly(torus300_incolle_a6[[k]][[1]], 2, 3, 10))
    return(append(aggr, list(time=time)))
    
  })
  save2Rdata(trs300_incole23_a6_aggrs)
}

###6点追加4,5セット目推定

{
  trs300_incole45_a6_aggrs<-lapply(4:5, function(k){
    
    cat("list", k, "calc\n")
    time<-system.time(aggr<-proposedMethodOnly(torus300_incolle_a6[[k]][[1]], 2, 3, 10))
    return(append(aggr, list(time=time)))
    
  })
  save2Rdata(trs300_incole45_a6_aggrs)
}


#補間前後のランドスケープ比較
##最初の1点を通る補間手法
#300点トーラスで比較

#補間前に不正解なのと補間後に正解したデータセット比較
in_crrct<-which(torus300_incolle13_aggrs[[1]][[2]] >= 0.5 & torus300_incolle13_aggrs[[1]][[2]] < 1.5)

wrng<-which(torus300_colle_aggrs[[1]][[2]] < 0.5 | torus300_colle_aggrs[[1]][[2]] >= 1.5)

in_crrct_wrng<-intersect(in_crrct, wrng)

##補間前不正解PD
trs300_1_w1_10_pd<-lapply(in_crrct_wrng[1:10], function(k)ripsDiag(torus300_colle_set[[1]][[k]][["noizyX"]], 2, 3, printProgress = T))
plotPDs(trs300_1_w1_10_pd)

##補間後正解PD
trs300_in1_w1_10_pd<-lapply(in_crrct_wrng[1:10], function(k)ripsDiag(torus300_incolle_set[[1]][[k]][["noizyX"]], 2, 3, printProgress = T))
plotPDs(trs300_in1_w1_10_pd)

##補間前不正解PL
trs300_1_w1_10_pls<-lapply(1:10, function(k)calc_landscape(trs300_1_w1_10_pd[[k]], 3))
par(cex.lab=2, cex.main=2, cex.axis=2, plt = c(0.2, 0.9, 0.2, 0.9), mfrow=c(2, 5), mgp=c(3.5,1,0), lwd=2)
plot_2ndpls(trs300_1_w1_10_pls, vert = F)

##補間後正解PL

trs300_in1_w1_10_pls<-lapply(1:10, function(k)calc_landscape(trs300_in1_w1_10_pd[[k]], 3))
plot_2ndpls(trs300_in1_w1_10_pls, vert = F)

#補間点の誤差を調べる
##最初の1点を通る補間手法
oldpar <- par(no.readonly=T)

in_trs300_1_w1_10errs<-lapply(in_crrct_wrng, function(i)torus_disterror(torus300_incolle_set[[1]][[i]][["noizyX"]], maxr = 2.5, minr = 1, nps = 300))
par(mgp=c(2.4,1,0))
boxplot(in_trs300_1_w1_10errs[1:10], xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6, lwd=2)




#時間と精度の関係グラフを作成
##1,2,3,4,6点追加を比較
##成功率をまとめる
a0_rate<-aggr_success_rates(torus300_colle_aggrs, c(2,1)) %>% do.call(rbind, .)

a1_rate<-aggr_success_rates(c(list(trs300_incolle1_a1_aggr2), trs300_incole13_a1_aggrs, trs300_incole45_a1_aggrs), c(2, 1)) %>% 
         do.call(rbind, .)

a2_rate<-aggr_success_rates(c(list(trs300_incolle1_a2_aggr), trs300_incole13_a2_aggrs, trs300_incole45_a2_aggrs), c(2, 1)) %>% 
  do.call(rbind, .)

a3_rate<-aggr_success_rates(c(list(trs300_incolle1_a3_aggr2), trs300_incole23_a3_aggrs, trs300_incole45_a3_aggrs), c(2, 1)) %>% 
  do.call(rbind, .)

a4_rate<-aggr_success_rates(c(list(trs300_incolle1_a4_aggr), trs300_incole23_a4_aggrs, trs300_incole45_a4_aggrs), c(2, 1)) %>% 
  do.call(rbind, .)

a6_rate<-aggr_success_rates(c(list(trs300_incolle1_a6_aggr), trs300_incole23_a6_aggrs, trs300_incole45_a6_aggrs), c(2, 1)) %>% 
  do.call(rbind, .)

all_rate<-append(aggr_success_rates(torus300_incolle13_aggrs, c(2, 1)), aggr_success_rates(torus300_incolle45_aggrs, c(2, 1))) %>% 
  do.call(rbind, .)


##計算時間をまとめる
##補間時間＋推定時間
a0_time<-sapply(torus300_colle_aggrs, function(t)t[["time"]][3])

a1_time<-c((torus300_incolle_a1[[1]][["time"]][3]+intt300_1_a1_time2)[3], sapply(torus300_incolle_a1, function(t)t[["time"]][3])[2:5]+sapply(c(trs300_incole13_a1_aggrs, trs300_incole45_a1_aggrs), function(t)t[["time"]][3]))

a2_time<-c((t300_a2_intime+intt300_1_a2_time)[3], sapply(torus300_incolle_a2, function(t)t[["time"]][3])[2:5]+sapply(c(trs300_incole13_a2_aggrs, trs300_incole45_a2_aggrs), function(t)t[["time"]][3]))

a3_time<-c((torus300_incolle_a3[[1]][["time"]][3]+intt300_1_a3_time2)[3], sapply(torus300_incolle_a3, function(t)t[["time"]][3])[2:5]+sapply(c(trs300_incole23_a3_aggrs, trs300_incole23_a3_aggrs),function(t)t[["time"]][3]))

a4_time<-c((t300_a4_intime+intt300_1_a4_time)[3], sapply(torus300_incolle_a4, function(t)t[["time"]][3])[2:5]+sapply(c(trs300_incole23_a4_aggrs, trs300_incole45_a4_aggrs), function(t)t[["time"]][3]))

a6_time<-c((t300_a6_intime+intt300_1_a6_time)[3], sapply(torus300_incolle_a6, function(t)t[["time"]][3])[2:5]+sapply(c(trs300_incole23_a6_aggrs, trs300_incole45_a6_aggrs), function(t)t[["time"]][3]))

all_time<-sapply(torus300_incolle_set_test, function(t)t[[2]][3])+c(sapply(torus300_incolle13_aggrs, function(t)t[["time"]][3]), sapply(torus300_incolle45_aggrs, function(t)t[["time"]][3]))

##図表作成
plot(c(a0_time, a1_time, a2_time, a3_time, a4_time, a6_time, all_time), c(unlist(a0_rate[,2]), unlist(a1_rate[,2]), unlist(a2_rate[,2]), unlist(a3_rate[,2]), unlist(a4_rate[,2]), unlist(a6_rate[,2]), unlist(all_rate[,2])), type="n",
     xlab="Computational time [sec]", ylab="Success Rates", cex.axis=1.6, cex.lab=1.6, xlim=c(0, 25000))
points(a0_time, unlist(a0_rate[,2]), col=1, pch=16, cex=1.5)
points(a1_time, unlist(a1_rate[,2]), col=2, pch=16, cex=1.5)
points(a2_time, unlist(a2_rate[,2]), col=3, pch=16, cex=1.5)
points(a3_time, unlist(a3_rate[,2]), col=4, pch=16, cex=1.5)
points(a4_time, unlist(a4_rate[,2]), col=5, pch=16, cex=1.5)
#points(a6_time, unlist(a6_rate[,2]), col="orange", pch=16, cex=1.5)
points(all_time, unlist(all_rate[,2]), col=6, pch=16, cex=1.5)
legend("bottomright", legend = c(sapply(0:4, function(k)paste0("P=", k)), "all vertexes"), col = c(1:6), pch = 16, cex=1.8)

rate_list<-list(a0_rate=a0_rate[,2],
                a1_rate=a1_rate[,2],
                a2_rate=a2_rate[,2],
                a3_rate=a3_rate[,2],
                a4_rate=a4_rate[,2],
                #a6_rate=a6_rate[,2],
                all_rate=all_rate[,2]
                )

time_list<-list(a0_time=a0_time,
                a1_time=a1_time,
                a2_time=a2_time,
                a3_time=a3_time,
                a4_time=a4_time,
                #a6_time=a6_time,
                all_time=all_time
                )

###点数増加倍率のラベルを付ける
library(maptools)
pointLabel(a0_time, unlist(a0_rate[,2]), as.character(rep(300, times=5)), cex=1.8)
pointLabel(a1_time, unlist(a1_rate[,2]), as.character(round(a1_points, digits = 2)), cex=1.8)
pointLabel(a2_time, unlist(a2_rate[,2]), as.character(round(a2_points, digits = 2)), cex=1.8)
pointLabel(a3_time, unlist(a3_rate[,2]), as.character(round(a3_points, digits = 2)), cex=1.8)
pointLabel(a4_time, unlist(a4_rate[,2]), as.character(round(a4_points, digits = 2)), cex=1.8)
#pointLabel(unlist(a6_rate[,2]), 1/a6_time,  as.character(round(a6_points[1:3], digits = 2)))
pointLabel(all_time, unlist(all_rate[,2]), as.character(round(all_points, digits = 2)), cex=1.8)

#全点に対して点数ラベルを付ける
pointLabel(c(a0_time, a1_time, a2_time, a3_time, a4_time, all_time), c(unlist(a0_rate[,2]), unlist(a1_rate[,2]), unlist(a2_rate[,2]), unlist(a3_rate[,2]), unlist(a4_rate[,2]), unlist(all_rate[,2])), c(as.character(rep(300, times=5)), 
           as.character(round(a1_points, digits = 2)), as.character(round(a2_points, digits = 2)), as.character(round(a3_points, digits = 2)), as.character(round(a4_points, digits = 2)), as.character(round(all_points, digits = 2))), cex=1.8)

###図表における平均・標準偏差計算
###平均
rate_mean<-sapply(rate_list, function(rate)mean(unlist(rate)))
time_mean<-sapply(time_list, function(time)mean(time))
lines(time_mean, rate_mean, lwd=2)

###標準偏差
rate_sd<-sapply(rate_list, function(rate)sd(unlist(rate)))
time_sd<-sapply(time_list, function(time)sd(time))
lines(rate_mean+rate_sd, 1/(time_mean+time_sd), lty="dashed")
lines(rate_mean-rate_sd, 1/(time_mean-time_sd), lty="dashed")

##点数増加率をまとめる
a1_points<-sapply(torus300_incolle_a1, function(trs){
  
  return(sapply(trs[[1]], function(X)X[["nsample"]]) %>%  mean())
  
}) 

a2_points<-sapply(torus300_incolle_a2, function(trs){
  
  return(sapply(trs[[1]], function(X)X[["nsample"]]) %>%  mean())
  
}) 

a3_points<-sapply(torus300_incolle_a3, function(trs){
  
  return(sapply(trs[[1]], function(X)X[["nsample"]]) %>%  mean())
  
}) 

a4_points<-sapply(torus300_incolle_a4, function(trs){
  
  return(sapply(trs[[1]], function(X)X[["nsample"]]) %>%  mean())
  
}) 

a6_points<-sapply(torus300_incolle_a6, function(trs){
  
  return(sapply(trs[[1]], function(X)X[["nsample"]]) %>% mean())
  
}) 

all_points<-sapply(torus300_incolle_set, function(trs){
  
  return(sapply(trs, function(X)X[["nsample"]]) %>% mean())
  
}) 

## plotly試し

torus350<-torus350_incolle_set[[2]][[26]][["noizyX"]] %>%
  dplyr::as_data_frame() %>% 
  cbind(., c(rep(1, 300), rep(2, 279)))
colnames(torus350)[4]<-"inter"
plotly::plot_ly(torus350, x = ~x, y = ~y, z = ~z, size = 1, color= ~inter, colors = c('#BF382A', '#0C4B8E')) %>% 
plotly::add_markers()
