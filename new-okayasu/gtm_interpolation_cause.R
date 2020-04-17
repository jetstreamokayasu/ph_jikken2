#GTMによる補間でなぜうまくいかないのか解明
##推定失敗データのランドスケープを見る

##補間点のみ点数削減に対してベッチ数推定
#300点トーラス、1~20セットのPD計算
trs300_incolle_1D_pd<-lapply(trs300_incolle_set1D[1:20], function(X)ripsDiag(X[["noizyX"]], 2, 3, printProgress = T))

#300点トーラス、1~10セットのベッチ数推定
{
  trs300_incolle_1D_1_10_aggr<-proposedMethodOnly(trs300_incolle_set1D[1:10], 2, 3, 10)
  save2Rdata(trs300_incolle_1D_1_10_aggr)
}

#300点トーラスのベッチ数推定
{
  trs300_incolle_1D_aggr<-proposedMethodOnly(trs300_incolle_set1D[1:10], 2, 3, 10)
  save2Rdata(trs300_incolle_1D_1_10_aggr)
}

#300点トーラス、1~20セットのパーシステントランドスケープ計算
trs300_incolle_1D_1_10_pls<-lapply(trs300_incolle_1D_pd, function(X)calcLandscape(X))
plot_lands(trs300_incolle_1D_1_10_pls[1:10], 1)
plot_lands(trs300_incolle_1D_1_10_pls[11:20], 1)

#300点トーラス、11~20セットのベッチ数推定
{
  trs300_incolle_1D_11_20_aggr<-proposedMethodOnly(trs300_incolle_set1D[11:20], 2, 3, 10)
  save2Rdata(trs300_incolle_1D_11_20_aggr)
}




##元の点と補間点両方を対象に点数削減
#300点トーラス、1~20セットのベッチ数推定
{
  trs300_incolle_1C_1_20_aggr<-proposedMethodOnly(trs300_incolle_set1C[1:20], 2, 3, 10)
  save2Rdata(trs300_incolle_1C_1_20_aggr)
}

#300点トーラス、1~20セットのPD計算
trs300_incolle_1C_pd<-lapply(trs300_incolle_set1C[1:20], function(X)ripsDiag(X[["noizyX"]], 2, 3, printProgress = T))

#300点トーラス、1~20セットのパーシステントランドスケープ計算
trs300_incolle_1C_1_20_pls<-lapply(trs300_incolle_1C_pd, function(X)calcLandscape(X))
plot_lands(trs300_incolle_1C_1_20_pls[1:10], 1)
plot_lands(trs300_incolle_1C_1_20_pls[11:20], 1)
