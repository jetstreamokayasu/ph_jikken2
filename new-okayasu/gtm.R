#GTM試し
library(myimg)
library(myfs)
library(TDA)
library(rgl)
library(phacm)
library(interpo3d)
library(seephacm)
require(phacm)
require(pracma)
require(deldir)
require(ggplot2)
require(plyr)
require(reshape2)
require(ggmap)
require(tidyverse)


##GTMによる次元削減
torus18_dist<-dist(torus.collect18[[1]][["noizyX"]])
trs18_vics1<-interpo3d:::get_vicinity(torus18_dist, 1, 30)
figurePlot3d(torus.collect18[[1]][["noizyX"]][-trs18_vics1, ])
points3d(torus.collect18[[1]][["noizyX"]][trs18_vics1, ], col=2)

trs18<-torus.collect18[[1]][["noizyX"]]

trs18_pca1<-prcomp(trs18[trs18_vics1, ])
plot(trs18_pca1[["x"]][, 1:2], col=4, pch=16)


library(maptools)
pointLabel(trs18_pca1[["x"]][, 1:2], as.character(trs18_vics1))

#GTM補間後、ベッチ数推定してみる
trs300_incolle_set1<-gtm_interpolate(torus300_colle_set[[1]][1:5], 30)
figurePlot3d(trs300_incolle_set1[[1]][["noizyX"]][1:300,])
points3d(trs300_incolle_set1[[1]][["noizyX"]][301:884,], col=2)

rgl.postscript("./data/trs300.eps", fmt="eps" )
rgl.postscript("./data/in_trs300.eps", fmt="eps" )


{
trs300_incolle_set1_test_aggr<-proposedMethodOnly(trs300_incolle_set1, 2, 3, 10)
save2Rdata(trs300_incolle_set1_test_aggr)
}


#GTM, 全ボロノイ領域の頂点補間, 点数削減

#300点トーラス
trs300_incolle_set1b<-gtm_inter_reduce(collect = torus300_colle_set[[1]], nvic = 30, ratio = 0.7)

{
  trs300_incolle_set1b_test_aggr<-proposedMethodOnly(trs300_incolle_set1b, 2, 3, 10)
  save2Rdata(trs300_incolle_set1b_test_aggr)
}

##300点トーラスで補間前、補間後、点数削減後を図表化
figurePlot3d(torus300_colle_set[[1]][[1]][["noizyX"]])
rgl.postscript("./data/trs300B.eps", fmt="eps")

trs330_1_1_inted<-voronoi_gtm_interpo(torus300_colle_set[[1]][[1]][["noizyX"]], nvics = 30)
trs300_1_1_upsam<-rbind(torus300_colle_set[[1]][[1]][["noizyX"]], trs300_1_1_inted)
points3d(trs300_1_1_inted, col=2)
rgl.postscript("./data/trs300B_inted.eps", fmt="eps")

trs300_1_1_red<-reduce_intered(intered_X = trs300_1_1_upsam, ratio = 0.95, n_ori = nrow(torus300_colle_set[[1]][[1]][["noizyX"]]))
figurePlot3d(trs300_1_1_red[["y"]][(trs300_1_1_red[["remain"]] <= 300), ])
rgl.postscript("./data/trs300B_ori.eps", fmt="eps")
points3d(trs300_1_1_red[["y"]][(trs300_1_1_red[["remain"]] > 300), ], col=2)
rgl.postscript("./data/trs300B_red.eps", fmt="eps")

#310点トーラス
trs310_incolle_set1<-gtm_inter_reduce(collect = torus310_colle_set[[1]], nvic = 30, ratio = 0.7)

{
  trs310_incolle_set1_test_aggr<-proposedMethodOnly(trs310_incolle_set1, 2, 3, 10)
  save2Rdata(trs310_incolle_set1_test_aggr)
}

trs310_incolle_set1<-gtm_inter_reduce(collect = torus310_colle_set[[1]], nvic = 30, ratio = 0.7)

#320点トーラス
trs320_incolle_set1<-gtm_inter_reduce(collect = torus320_colle_set[[1]], nvic = 30, ratio = 0.7)
{
  trs320_incolle_set1_test_aggr<-proposedMethodOnly(trs320_incolle_set1, 2, 3, 10)
  save2Rdata(trs320_incolle_set1_test_aggr)
}

#330点トーラス
trs330_incolle_set1<-gtm_inter_reduce(collect = torus330_colle_set[[1]], nvic = 30, ratio = 0.7)
{
  trs330_incolle_set1_test_aggr<-proposedMethodOnly(trs330_incolle_set1, 2, 3, 10)
  save2Rdata(trs330_incolle_set1_test_aggr)
}

{
#340点トーラス
trs340_incolle_set1<-gtm_inter_reduce(collect = torus340_colle_set[[1]], nvic = 30, ratio = 0.95)
{
  trs340_incolle_set1_test_aggr<-proposedMethodOnly(trs340_incolle_set1, 2, 3, 10)
  save2Rdata(trs340_incolle_set1_test_aggr)
}
}

#350点トーラス

trs350_incolle_set1<-gtm_inter_reduce(collect = torus350_colle_set[[1]], nvic = 30, ratio = 0.7)
{
  trs350_incolle_set1_test_aggr<-proposedMethodOnly(trs350_incolle_set1, 2, 3, 10)
  save2Rdata(trs350_incolle_set1_test_aggr)
}

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

##並列化を試す
##350点トーラス
int_time<-system.time(trs350_incolle_set1b<-parLapply(cl, torus350_colle_set[[1]], function(X){
  
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
##350点トーラス
{
  trs350_incolle_set1b_test_aggr<-proposedMethodOnly(trs350_incolle_set1b, 2, 3, 10)
  save2Rdata(trs350_incolle_set1b_test_aggr)
}


##並列化を試す
##300点トーラス


int_time300p<-system.time(trs300_incolle_set1c<-parLapply(cl, torus300_colle_set[[1]], function(X){
  
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

#parallelで補間
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


##不具合チェック
trs340_set89_inted<-voronoi_gtm_interpo(torus340_colle_set[[1]][[89]][["noizyX"]], nvics = 30)

trs340_1_89_vic10s<-interpo3d:::get_vicinity(dist(torus340_colle_set[[1]][[89]][["noizyX"]]), 10, 30)

#deldirエラーチェック
##Both vertex orderings are clockwise. See help for deldir.
trs350_1_89_inted<-voronoi_gtm_interpo(torus350_colle_set[[1]][[89]][["noizyX"]], nvics = 30)

trs350_1_89_vic39s<-interpo3d:::get_vicinity(dist(torus350_colle_set[[1]][[89]][["noizyX"]]), 39, 30)

X<-torus340_colle_set[[1]][[89]][["noizyX"]][trs340_1_89_vic10s,]
# (分散が0の変数を削除した後に) 1. オートスケーリング
Var0Variable <- which(apply(X,2,var) == 0)
if (length(Var0Variable) == 0) {
  #print("分散が0の変数はありません")
} else {
  sprintf("分散が0の変数が %d つありました", length(Var0Variable))
  print( "変数:" )
  print( Var0Variable )
  print( "これらを削除します" )
  X <- X[,-Var0Variable]
}
X <- scale(X, center = TRUE, scale = TRUE)
# 2. マップサイズ
MapsizeColumn = 15 #横 10
MapsizeRow = 15 #縦 10 、
# 3. 動径基底関数 (Radial Basis Function, RBF) の数
RBFsizeColumn = 3 #横 3
RBFsizeRow = 3 #縦 3
# 4. ガウス関数の分散
RBFVariance = 1
# 5. EMアル
# 6. データ空間における分散の逆数βの初期値
# 7. 重みWの初期値
XGrid = gtm.rctg( MapsizeColumn, MapsizeRow)#写像先のグリッド
RBFGrid = gtm.rctg( RBFsizeColumn, RBFsizeRow)#基底関数のグリッド
RBFSetup = gtm.gbf( RBFGrid, RBFVariance^(1/2), XGrid)
InitnalWBeta = gtm.pci.beta(X, XGrid, RBFSetup)
#Beta = 0.01
Beta = InitnalWBeta$beta
# 9. GTMマップ作成
GTMResults = gtm.trn( X, RBFSetup, InitnalWBeta$W, Lambda, NumOfTraining, Beta, quiet = T)
# 10. 二次元のマップ上でサンプルの位置関係を確認
GTMDist = gtm.dist(X, RBFSetup %*% GTMResults$W)
GTMR = gtm.resp3(GTMDist, GTMResults$beta, ncol(X))$R
GTMMean = t(GTMR) %*% XGrid

par(pty="s")
plot(GTMMean[,1], GTMMean[,2])
points(GTMMean[1,1], GTMMean[1,2], pch=16, col=2)
#次元削減後のボロノイ分割
res<-try(deldir(GTMMean[,1], GTMMean[,2], eps = 1e-5))
trs340_1_89_tiles158<-tile.list(res)
for(i in 1:res$n.data){	polygon(trs340_1_89_tiles158[[i]], lwd=2) }

#境界領域に接しないボロノイ領域の頂点を選ぶ
vertx<-cbind(res[["dirsgs"]][["x1"]][!res[["dirsgs"]][["bp1"]]], res[["dirsgs"]][["y1"]][!res[["dirsgs"]][["bp1"]]])
points(vertx, pch=16, col=4)
RBF_inter = gtm.gbf(RBFGrid, RBFVariance^(1/2), vertx)
#inter_dist<-gtm.dist(X, RBF_inter %*% GTMResults$W)
# inter_R = gtm.resp3(inter_dist, GTMResults$beta, ncol(X))$R
# inter_mean = t(inter_R) %*% (RBF_inter %*% GTMResults$W)
inter_mean<-RBF_inter %*% GTMResults$W
inter_inv<-sapply(1:nrow(inter_mean), function(k){inter_mean[k, ] * attr(X, "scaled:scale") + attr(X, "scaled:center")})
inter_inv<-t(inter_inv)

#補間後点数削減チェック
trs340_1_89_upsam<-rbind(torus340_colle_set[[1]][[89]][["noizyX"]], trs340_set89_inted)

trs340_1_89_red<-reduce_intered(intered_X = trs340_1_89_upsam, ratio = 0.95, n_ori = nrow(torus340_colle_set[[1]][[89]][["noizyX"]]))

thresh<-quantile_threshold(x = 0.95, X = trs340_1_89_upsam)
trs340_1_89_cell<-cell_set2(x = trs340_1_89_upsam, thresh = thresh)



