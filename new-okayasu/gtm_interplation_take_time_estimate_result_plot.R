#並列化準備
#parallel使用
library(parallel)

cl <- makeCluster(4, outfile="")

clusterEvalQ(cl,{
  library(phacm)
  library(interpo3d)
  library(pracma)
  library(deldir)
  library(gtm)
  library(tidyverse)
  library(myfs)
  library(TDA)
})

clusterEvalQ(cl, {
  source('~/R/interpolation_test/interpo_func.R', encoding = 'UTF-8')
  source('~/R/interpolation_test/interpo_func_2.R', encoding = 'UTF-8')
  source('~/R/interpolation_test/interpo_func.R', encoding = 'UTF-8')
  source('~/R/p_reduce/reduce_func.R', encoding = 'UTF-8')
  source('~/R/interpolation_test/reduce_func2.R', encoding = 'UTF-8')
}
)


stopCluster(cl)


#補間時間を測定
#元の点と補間点両方対象の点数削減、補間点のみの点数削減でベッチ数推定制度比較

#元のデータ点はすべて残し、補間された点のみ減らす
#300点トーラス

int_time300p2<-system.time(trs300_incolle_set1D<-parLapply(cl, torus300_colle_set[[1]], function(X){
  
  inter_oricord<-voronoi_gtm_interpo(X[[2]], nvics = 30)
  inter_oricord<-inter_oricord[!is.na(inter_oricord[,1]), ]
  red_oricord<-reduce_intered2(intered_X = inter_oricord, ratio = 0.99)
  X[[2]]<-rbind(X[[2]], red_oricord)
  X[[1]]<-nrow(X[[2]])
  sink(paste0("./parallel/", X[[1]], "data", format(Sys.time(), "%m%d_%H%M"), ".txt"))
  print(paste("dataset has", X[[1]], "points"))
  sink()
  return(X)
  
}))

#parallelで補間後ベッチ数推定
##300点トーラス
{
  trs300_incolle_set1D_test_aggr<-proposedMethodOnly(trs300_incolle_set1D, 2, 3, 10)
  save2Rdata(trs300_incolle_set1D_test_aggr)
}



#元の点および補間点を両方を対象として点数削減
##300点トーラス
int_time300p<-system.time(trs300_incolle_set1C<-parLapply(cl, torus300_colle_set[[1]], function(X){
  
  inter_oricord<-voronoi_gtm_interpo(X[[2]], nvics = 30)
  inter_oricord<-inter_oricord[!is.na(inter_oricord[,1]), ]
  red_oricord<-reduce_intered(intered_X = rbind(X[[2]], inter_oricord), ratio = 0.95, n_ori = nrow(X[[2]]))
  X[[2]]<-red_oricord[["y"]]
  X[[1]]<-nrow(X[[2]])
  sink(paste0("./parallel/", X[[1]], "data", format(Sys.time(), "%m%d_%H%M"), ".txt"))
  print(paste("dataset has", X[[1]], "points"))
  sink()
  return(X)
  
}))

#parallelで補間後ベッチ数推定
##300点トーラス
{
  trs300_incolle_set1C_test_aggr<-proposedMethodOnly(trs300_incolle_set1c, 2, 3, 10)
  save2Rdata(trs300_incolle_set1C_test_aggr)
}

