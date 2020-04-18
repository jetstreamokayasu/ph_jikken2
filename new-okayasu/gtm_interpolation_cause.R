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
  trs300_incolle_1D_aggr<-proposedMethodOnly(trs300_incolle_set1D, 2, 3, 10)
  save2Rdata(trs300_incolle_1D_aggr)
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

#300点トーラスのベッチ数推定
{
  trs300_incolle_1C_aggr<-proposedMethodOnly(trs300_incolle_set1C, 2, 3, 10)
  save2Rdata(trs300_incolle_1C_aggr)
}

#補間前の300点トーラスの1セット目のベッチ数推定
torus300_colset1_aggrs<-proposedMethodOnly(torus300_colle_set[[1]], 2, 3, 10)

#補間前の推定で正しいベッチ数を推定したトーラス
trs300_1_rdx<-range_index(torus300_colset1_aggrs[[1]], min = 1.5, max = 2.5)
  
#元の点と補間点を対象とした点数削減後に不正確に推定されたトーラス
trs300_incole_1C_widx<-range_index(trs300_incolle_1C_aggr[[1]], min = 0, max = 1.5)

#補間点のみを対象とした点数削減後に不正確に推定されたトーラス
trs300_incole_1D_widx<-range_index(trs300_incolle_1D_aggr[[1]], min = 0, max = 1.5)

#補間前に正しく、補間後に誤って推定されたトーラスのインデックス
#元の点と補間点を対象とした点数削減後
trs300_1C_intsec_idx<-intersect(trs300_1_rdx, trs300_incole_1C_widx)

#補間前に正しく、補間後に誤って推定されたトーラスのインデックス
#補間点のみを対象とした点数削減後
trs300_1D_intsec_idx<-intersect(trs300_1_rdx, trs300_incole_1D_widx)

#補間前に正しく、補間後に誤って推定されたトーラスのPD
#補間点のみを対象とした点数削減後
trs300_incolle_1D_intsec_pd<-parLapply(trs300_incolle_set1D[trs300_1D_intsec_idx], function(X){
  
  pd<-ripsDiag(X[["noizyX"]], 2, 3, printProgress = T)
  sink(paste0("./parallel/", X[[1]], "data", format(Sys.time(), "%m%d_%H%M"), ".txt"))
  print(paste("dataset has", X[[1]], "points"))
  sink()
  
  return(pd)
  
  })
