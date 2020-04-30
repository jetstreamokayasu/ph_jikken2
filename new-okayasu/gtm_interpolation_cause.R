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
trs300_incolle_1D_intsec_pd<-parLapply(cl, trs300_incolle_set1D[trs300_1D_intsec_idx], function(X){
  
  pd<-ripsDiag(X[["noizyX"]], 2, 3, printProgress = T)
  sink(paste0("./parallel/", X[[1]], "data", format(Sys.time(), "%m%d_%H%M"), ".txt"))
  print(paste("dataset has", X[[1]], "points"))
  sink()
  
  return(pd)
  
  })

#補間前に正しく、補間後に誤って推定されたトーラスのPD
#元の点と補間点を対象とした点数削減後
trs300_incolle_1C_intsec_pd<-lapply(trs300_incolle_set1C[trs300_1C_intsec_idx], function(X)ripsDiag(X[["noizyX"]], 2, 3, printProgress = T))

#補間前の300点トーラスのPD
trs300_colle1_pd<-lapply(torus300_colle_set[[1]], function(X)ripsDiag(X[["noizyX"]], 2, 3, printProgress = T))


#補間前に正しく、全点対象、補間点のみそれぞれの点数削減の時に誤って推定されたインデックス
trs300_1C_1D_intsec_idx<-intersect(trs300_1D_intsec_idx, trs300_1C_intsec_idx)

#補間前に正しく、全点対象および補間点のみそれぞれの点数削減の時に誤って推定されたデータセットのPL
#補間前10セット
trs300_colle1_pls<-lapply(trs300_colle1_pd[trs300_1C_1D_intsec_idx[1:10]], function(X)calcLandscape(X))
plot_lands(trs300_colle1_pls, 1)

#補間後10セット。全点対象点数削減後
#補間前正解、補間＋全点対象点数削減後不正解のインデックス
trs300_1C_intsec_idx2<-trs300_1C_intsec_idx %in% trs300_1C_1D_intsec_idx[1:10]

#補間前正解、補間＋全点対象点数削減後不正解のPL
trs300_incolle_1C_intsec_pls<-lapply(trs300_incolle_1C_intsec_pd[trs300_1C_intsec_idx2], function(X)calcLandscape(X))
plot_lands(trs300_incolle_1C_intsec_pls, 1)

#補間後10セット。補間点のみ対象点数削減後
#補間前正解、補間＋補間点のみ対象点数削減後不正解のインデックス
trs300_1D_intsec_idx2<-trs300_1D_intsec_idx %in% trs300_1C_1D_intsec_idx[1:10]

#補間前正解、補間＋補間点のみ対象点数削減後不正解のPL
trs300_incolle_1D_intsec_pls<-lapply(trs300_incolle_1D_intsec_pd[trs300_1D_intsec_idx2], function(X)calcLandscape(X))
plot_lands(trs300_incolle_1D_intsec_pls, 1)

#誤って推定されたデータセットの補間点の誤差を調べる
#補間点のみ削除の場合
trs300_incolle_1D_errs<-lapply(trs300_1C_1D_intsec_idx[1:10], function(i)torus_disterror(trs300_incolle_set1D[[i]][["noizyX"]], maxr = 2.5, minr = 1, nps = 300))
oldpar <- par(no.readonly=T)
par(mgp=c(2.4,1,0))
boxplot(trs300_incolle_1D_errs, xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6, lwd=2)

#全点対象点数削減
trs300_incolle_1C_errs<-lapply(trs300_1C_1D_intsec_idx[1:10], function(i)torus_disterror(trs300_incolle_set1C[[i]][["noizyX"]], maxr = 2.5, minr = 1, nps = 0))
oldpar <- par(no.readonly=T)
par(mgp=c(2.4,1,0))
boxplot(trs300_incolle_1C_errs, xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6, lwd=2)


#---------------------------------------

#補間前後でどちらも正しく推定されたデータセットのPL比較

#元の点と補間点を対象とした点数削減後に正確に推定されたトーラス
trs300_incole_1C_ridx<-range_index(trs300_incolle_1C_aggr[[1]], min = 1.5, max = 2.5)

#補間点のみを対象とした点数削減後に正確に推定されたトーラス
trs300_incole_1D_ridx<-range_index(trs300_incolle_1D_aggr[[1]], min = 1.5, max = 2.5)

#補間前後でどちらも正しく推定されたトーラスのインデックス
trs300_1C_1D_intsec_ridx<-intersect(trs300_1_rdx, trs300_incole_1C_ridx) %>% intersect(., trs300_incole_1D_ridx)

#補間前後でどちらも正しく推定されたトーラスのPD
#全点対象点数削減
trs300_incolle_1C_ridx_pd<-lapply(trs300_incolle_set1C[trs300_1C_1D_intsec_ridx], function(X)ripsDiag(X[["noizyX"]], 2, 3, printProgress = T))

#補間前後でどちらも正しく推定されたトーラスのPD
#補間点のみ点数削減
trs300_incolle_1D_ridx_pd<-lapply(trs300_incolle_set1D[trs300_1C_1D_intsec_ridx], function(X)ripsDiag(X[["noizyX"]], 2, 3, printProgress = T))

#補間前後でどちらも正しく推定されたトーラスのPL
#補間前10セット
trs300_colle1_pls_B<-lapply(trs300_colle1_pd[trs300_1C_1D_intsec_ridx[1:10]], function(X)calc_landscape(X, maxscale = 3))
plot_lands(trs300_colle1_pls_B, 1)

#全点対象点数削減
trs300_incolle_1C_ridx_pls<-lapply(trs300_incolle_1C_ridx_pd, function(X)calc_landscape(X, maxscale = 3))
plot_lands(trs300_incolle_1C_ridx_pls[1:10], 1)

#補間点のみ削減
trs300_incolle_1D_ridx_pls<-lapply(trs300_incolle_1D_ridx_pd, function(X)calc_landscape(X, maxscale = 3))
plot_lands(trs300_incolle_1D_ridx_pls[1:10], 1)

#正しく推定されたデータセットの補間点の誤差を調べる
#補間点のみ点数削減
trs300_incolle_1D_errs_B<-lapply(trs300_1C_1D_intsec_ridx[1:10], function(i)torus_disterror(trs300_incolle_set1D[[i]][["noizyX"]], maxr = 2.5, minr = 1, nps = 300))
par(mgp=c(2.4,1,0))
boxplot(trs300_incolle_1D_errs_B, xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6, lwd=2)

#全点対象点数削減
trs300_incolle_1C_errs_B<-lapply(trs300_1C_1D_intsec_ridx[1:10], function(i)torus_disterror(trs300_incolle_set1C[[i]][["noizyX"]], maxr = 2.5, minr = 1, nps = 0))
boxplot(trs300_incolle_1C_errs_B, xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6, lwd=2)

