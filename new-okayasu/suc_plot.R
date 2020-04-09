#成功率のグラフ化

#プロット領域準備
plot(sucs_rates_tidy, pch=16, cex.axis=1.6, xlab="Data Density", ylab="Success Rates", cex.lab=1.6, xlim=c(300, 350), ylim=c(0.0, 1.0), xaxt="n", type="n")
axis(side=1, at=seq(300, 350, by=10), labels=c(paste0(seq(30, 35), "/(pi^2)")), cex.axis=1.1)

points(rep(300, 5), sucs300_rate[,2], pch=16, col='#7F7F7F')
points(rep(310, 5), sucs310_rate[,2], pch=16, col='#7F7F7F')
points(rep(320, 5), sucs320_rate[,2], pch=16, col='#7F7F7F')
points(rep(330, 5), sucs330_rate[,2], pch=16, col='#7F7F7F')
points(rep(340, 5), sucs340_rate[,2], pch=16, col='#7F7F7F')
points(rep(350, 5), sucs350_rate[,2], pch=16, col='#7F7F7F')

sucsrate_mean<-sapply(sucs_rates, function(rate)mean(unlist(rate[,2])))
lines(seq(300, 350, by=10), sucsrate_mean, col='#7F7F7F')

sucsrate_sd<-sapply(sucs_rates, function(rate)sd(unlist(rate[,2])))
lines(seq(300, 350, by=10), sucsrate_mean-sucsrate_sd, lty="dashed", col='#7F7F7F')
lines(seq(300, 350, by=10), sucsrate_mean+sucsrate_sd, lty="dashed", col='#7F7F7F')

##補間前1次ベッチ数
plot(sucs_rates_tidy, pch=16, cex.axis=1.6, xlab="Data Density", ylab="Success Rates", cex.lab=1.6, xlim=c(300, 350), ylim=c(0.0, 1.0), xaxt="n", type="n")
axis(side=1, at=seq(300, 350, by=10), labels=c(paste0(seq(30, 35), "/(pi^2)")), cex.axis=1.1)

points(rep(300, 5), sucs300_rate[,1], pch=16, col='#7F7F7F')
points(rep(310, 5), sucs310_rate[,1], pch=16, col='#7F7F7F')
points(rep(320, 5), sucs320_rate[,1], pch=16, col='#7F7F7F')
points(rep(330, 5), sucs330_rate[,1], pch=16, col='#7F7F7F')
points(rep(340, 5), sucs340_rate[,1], pch=16, col='#7F7F7F')
points(rep(350, 5), sucs350_rate[,1], pch=16, col='#7F7F7F')

sucsrate_mean<-sapply(sucs_rates, function(rate)mean(unlist(rate[,1])))
lines(seq(300, 350, by=10), sucsrate_mean, col='#7F7F7F')

sucsrate_sd<-sapply(sucs_rates, function(rate)sd(unlist(rate[,1])))
lines(seq(300, 350, by=10), sucsrate_mean-sucsrate_sd, lty="dashed", col='#7F7F7F')
lines(seq(300, 350, by=10), sucsrate_mean+sucsrate_sd, lty="dashed", col='#7F7F7F')


#補間後
torus350_incolle_rate<-aggr_success_rates(torus350_incolle_aggrs, c(2,1))
torus340_incolle13_rate<-aggr_success_rates(torus340_incolle13_aggrs, c(2,1))
torus330_incolle13_rate<-aggr_success_rates(torus330_incolle13_aggrs, c(2,1))
torus320_incolle13_rate<-aggr_success_rates(torus320_incolle13_aggrs, c(2,1))
torus310_incolle13_rate<-aggr_success_rates(torus310_incolle13_aggrs, c(2,1))
torus300_incolle13_rate<-aggr_success_rates(torus320_incolle13_aggrs, c(2,1))

torus340_incolle45_rate<-aggr_success_rates(torus340_incolle45_aggrs, c(2,1))
torus330_incolle45_rate<-aggr_success_rates(torus330_incolle45_aggrs, c(2,1))
torus320_incolle45_rate<-aggr_success_rates(torus320_incolle45_aggrs, c(2,1))
torus310_incolle45_rate<-aggr_success_rates(torus320_incolle45_aggrs, c(2,1))
torus300_incolle45_rate<-aggr_success_rates(torus320_incolle45_aggrs, c(2,1))


in350_rates<-do.call(rbind, torus350_incolle_rate)
in340_rates<-do.call(rbind, append(torus340_incolle13_rate, torus340_incolle45_rate))
in330_rates<-do.call(rbind, append(torus330_incolle13_rate, torus330_incolle45_rate))
in320_rates<-do.call(rbind, append(torus320_incolle13_rate, torus320_incolle45_rate))
in310_rates<-do.call(rbind, append(torus310_incolle13_rate, torus310_incolle45_rate))
in300_rates<-do.call(rbind, append(torus300_incolle13_rate, torus300_incolle45_rate))


in300_1_rate<-do.call(rbind, torus300_1_1_rate)

insub_rates<-list("300"=in300_rates,
                  "310"=in310_rates,
                  "320"=in320_rates,
                  "330"=in330_rates,
                  "340"=in340_rates,
                  "350"=in350_rates)

#2次ベッチ数新手法補間後
points(rep(300, 5), in300_rates[,2], col='#4C4CFF', pch=16)
points(rep(310, 5), in310_rates[,2], col='#4C4CFF', pch=16)
points(rep(320, 5), in320_rates[,2], col='#4C4CFF', pch=16)
points(rep(330, 5), in330_rates[,2], col='#4C4CFF', pch=16)
points(rep(340, 5), in340_rates[,2], col='#4C4CFF', pch=16)
points(rep(350, 5), in350_rates[,2], col='#4C4CFF', pch=16)

in_dim2_mean<-sapply(insub_rates, function(rate)mean(unlist(rate[,2])))
lines(seq(300, 350, by=10), in_dim2_mean, col='#4C4CFF')

in_dim2_sd<-sapply(insub_rates, function(rate)sd(unlist(rate[,2])))
lines(seq(300, 350, by=10), in_dim2_mean-in_dim2_sd, lty="dashed", col='#4C4CFF')
lines(seq(300, 350, by=10), in_dim2_mean+in_dim2_sd, lty="dashed", col='#4C4CFF')


##補間後1次ベッチ数
points(rep(300, 5), in300_rates[,1], col='#4C4CFF', pch=16)
points(rep(310, 5), in310_rates[,1], col='#4C4CFF', pch=16)
points(rep(320, 5), in320_rates[,1], col='#4C4CFF', pch=16)
points(rep(330, 5), in330_rates[,1], col='#4C4CFF', pch=16)
points(rep(340, 5), in340_rates[,1], col='#4C4CFF', pch=16)
points(rep(350, 5), in350_rates[,1], col='#4C4CFF', pch=16)

in_dim2_mean<-sapply(insub_rates, function(rate)mean(unlist(rate[,1])))
lines(seq(300, 350, by=10), in_dim2_mean, col='#4C4CFF')

in_dim2_sd<-sapply(insub_rates, function(rate)sd(unlist(rate[,1])))
lines(seq(300, 350, by=10), in_dim2_mean-in_dim2_sd, lty="dashed", col='#4C4CFF')
lines(seq(300, 350, by=10), in_dim2_mean+in_dim2_sd, lty="dashed", col='#4C4CFF')


#スライド用トーラス画像

figurePlot3d(torus300_colle_set[[1]][[1]][["noizyX"]])
rgl.snapshot("./data/torus_img/torus0.png") 

for(j in 1:4){
  
  figurePlot3d(torus320_colle_set[[1]][[j]][["noizyX"]])
  rgl.snapshot(paste0("./data/torus_img/torus320_", j, ".png"))
  
}


#GTM補間および点数削減後の精度
##2次ベッチ数成功率
gtm_300rate_dim2<-cycle_number(trs300_incolle_set1b_test_aggr, 2)[2]/100
gtm_310rate_dim2<-cycle_number(trs310_incolle_set1_test_aggr, 2)[2]/100
gtm_320rate_dim2<-cycle_number(trs320_incolle_set1_test_aggr, 2)[2]/100
gtm_330rate_dim2<-cycle_number(trs330_incolle_set1_test_aggr, 2)[2]/100
gtm_340rate_dim2<-cycle_number(trs340_incolle_set1_test_aggr, 2)[2]/100

##1次ベッチ数成功率
gtm_300rate_dim1<-cycle_number(trs300_incolle_set1b_test_aggr, 1)[3]/100
gtm_310rate_dim1<-cycle_number(trs310_incolle_set1_test_aggr, 1)[3]/100
gtm_320rate_dim1<-cycle_number(trs320_incolle_set1_test_aggr, 1)[3]/100
gtm_330rate_dim1<-cycle_number(trs330_incolle_set1_test_aggr, 1)[3]/100
gtm_340rate_dim1<-cycle_number(trs340_incolle_set1_test_aggr, 1)[3]/100

#2次ベッチ数成功率をプロット
points(seq(300, 340, by=10), c(gtm_300rate_dim2, gtm_310rate_dim2, gtm_320rate_dim2, gtm_330rate_dim2, gtm_340rate_dim2), col=2, pch=16, cex=1.5)
lines(seq(300, 340, by=10), c(gtm_300rate_dim2, gtm_310rate_dim2, gtm_320rate_dim2, gtm_330rate_dim2, gtm_340rate_dim2), col=2, lwd=2)

#1次ベッチ数成功率をプロット
points(seq(300, 340, by=10), c(gtm_300rate_dim1, gtm_310rate_dim1, gtm_320rate_dim1, gtm_330rate_dim1, gtm_340rate_dim1), col=2, pch=16, cex=1.5)
lines(seq(300, 340, by=10), c(gtm_300rate_dim1, gtm_310rate_dim1, gtm_320rate_dim1, gtm_330rate_dim1, gtm_340rate_dim1), col=2, lwd=2)
