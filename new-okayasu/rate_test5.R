require(phacm)
require(tidyverse)
require(TDA)
require(myfs)
require(rgl)
library(seephacm)
library(tictoc)
library(devtools)
library(interpo3d)

#ボロノイ領域の中心点から最も遠い頂点に補うように変更
#成功率を調べる

intrs300_1<-all_interpolate2(collect = torus300_colle_set[[1]], nvic = 15)

intt300_1_time2<-system.time(intrs300_1_aggr2<-proposedMethodOnly(intrs300_1, 2, 3, 10))


#平均よりも遠い点に補間
intrs300_1_2nd<-all_interpolate2(collect = torus300_colle_set[[1]], nvic = 15)

intt300_1_time3<-system.time(intrs300_1_aggr3<-proposedMethodOnly(intrs300_1_2nd, 2, 3, 10))

#補間後の平均増加率を調べる
trs300_incolle1_points<-unlist(lapply(torus300_incolle_set[[1]], function(X)X[[1]])) %>% '/' (., 300) %>% mean() 
                        

#平均よりも遠い点に補間手法で成功率を調べる
#300点補間
t300_intime2<-system.time(torus300_incolle_set2<-lapply(torus300_colle_set, function(k)all_interpolate2(k, 15)))

##300点トーラス補間後1~3セット目を推定
torus300_incolle13_aggrs2<-lapply(1:3, function(k){
  
  cat("list", k, "calc\n")
  time<-system.time(aggr<-proposedMethodOnly(torus300_incolle_set2[[k]], 2, 3, 10))
  save(aggr, file = paste0("./data/in300_aggr2_", k, ".RData"))
  return(append(aggr, list(time=time)))
  
})
save2Rdata(torus300_incolle13_aggrs2)

#310点補間
t310_intime2<-system.time(torus310_incolle_set2<-lapply(torus310_colle_set, function(k)all_interpolate2(k, 15)))

##310点トーラス補間後1~3セット目を推定
torus310_incolle13_aggrs2<-lapply(1:3, function(k){
  
  cat("list", k, "calc\n")
  time<-system.time(aggr<-proposedMethodOnly(torus310_incolle_set2[[k]], 2, 3, 10))
  save(aggr, file = paste0("./data/in310_aggr2_", k, ".RData"))
  return(append(aggr, list(time=time)))
  
})
save2Rdata(torus310_incolle13_aggrs2)

##320点トーラス補間後1~3セット目を推定
torus320_incolle13_aggrs2<-lapply(1:3, function(k){
  
  cat("list", k, "calc\n")
  time<-system.time(aggr<-proposedMethodOnly(torus320_incolle_set[[k]], 2, 3, 10))
  save(aggr, file = paste0("./data/in310_aggr2_", k, ".RData"))
  return(append(aggr, list(time=time)))
  
})
save2Rdata(torus310_incolle13_aggrs2)