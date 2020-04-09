#補間後に点数を減らしてみる

thresh<-quantile_threshold(x = 0.2, X = torus300_incolle_set[[1]][[1]][["noizyX"]])

cell<-cell_set2(torus300_incolle_set[[1]][[1]][["noizyX"]], thresh = thresh)
cnct<-connect2(1, cell_p = cell, all = 1:nrow(torus300_incolle_set[[1]][[1]][["noizyX"]]))
trs_300_red<-reduce_points(x = torus300_incolle_set[[1]][[1]][["noizyX"]], conect = cnct)

trs_300_redc<-quantile_threshold(x = 0.2, X = torus300_incolle_set[[1]][[2]][["noizyX"]]) %>% 
              cell_set2(torus300_incolle_set[[1]][[2]][["noizyX"]], thresh = .) %>% 
              connect2(1, cell_p = ., all = 1:nrow(torus300_incolle_set[[1]][[2]][["noizyX"]])) %>% 
              reduce_points(x = torus300_incolle_set[[1]][[2]][["noizyX"]], conect = .)

# 減らした後のトーラスを描画。元のデータセットにない点を赤くする
figurePlot3d(torus300_colle_set[[1]][[1]][["noizyX"]])
points3d(trs_300_redc[which(!(trs_300_redc[,1] %in% torus300_colle_set[[1]][[2]][["noizyX"]][,1])), ], col=2)

# 300点トーラス100個の1セット目を試しに推定してみる
## 補間後の1セット目全体を各点間の距離の20パーセンタイルを閾値として点を減らす
trs300_incolle1_redc<-lapply(torus300_incolle_set[[1]], function(X){
  
  redc<-quantile_threshold(x = 0.2, X = X[["noizyX"]]) %>% 
    cell_set2(X[["noizyX"]], thresh = .) %>% 
    connect2(1, cell_p = ., all = 1:nrow(X[["noizyX"]])) %>% 
    reduce_points(x = X[["noizyX"]], conect = .)
  
  return(list(nsample=nrow(redc), noizyX=redc, diag=0))
  
})

trs300_1_redc_time<-system.time(trs300_incolle1_redc_aggr<-proposedMethodOnly(trs300_incolle1_redc, 2, 3, 10))

#ボロノイ領域の頂点の一部だけに補間した場合のデータ店の増加率を調べる
## 各データセットのデータ点数を抽出
trs300_incolle1_a1_points<-unlist(lapply(torus300_incolle_1_a1, function(X)X[[1]]))

trs300_incolle1_a3_points<-unlist(lapply(torus300_incolle_1_a3, function(X)X[[1]]))

trs300_incolle1_a6_points<-unlist(lapply(torus300_incolle_1_a6, function(X)X[[1]]))

trs300_incolle1_points<-unlist(lapply(torus300_incolle_set[[1]], function(X)X[[1]]))

trs300_incolle1_redc_points<-unlist(lapply(trs300_incolle1_redc, function(X)X[[1]]))
trs300_incolle1_redc_incre<-trs300_incolle1_redc_points/300
mean(trs300_incolle1_redc_incre-1)

## 増加率の平均を計算
trs300_incolle1_a1_incre<-trs300_incolle1_a1_points/300
mean(trs300_incolle1_a1_incre-1)

trs300_incolle1_a3_incre<-trs300_incolle1_a3_points/300
mean(trs300_incolle1_a3_incre-1)

trs300_incolle1_a6_incre<-trs300_incolle1_a6_points/300
mean(trs300_incolle1_a6_incre-1)

trs300_incolle1_incre<-trs300_incolle1_points/300
mean(trs300_incolle1_incre-1)


trs300_incolle1_a2_incre<-unlist(lapply(torus300_incolle_1_a2, function(X)X[[1]])) %>% 
                          '/'(., 300)
mean(trs300_incolle1_a2_incre-1)


trs300_incolle1_a4_incre<-unlist(lapply(torus300_incolle_1_a4, function(X)X[[1]])) %>% 
  '/'(., 300)
mean(trs300_incolle1_a4_incre-1)



#補間された点のみを減らす
#元のデータはすべて残す
trs350_1_1_inted<-voronoi_gtm_interpo(torus350_colle_set[[1]][[1]][["noizyX"]], nvics = 30)
trs350_1_1_inted<-trs350_1_1_inted[!is.na(trs350_1_1_inted[,1]), ]
trs350_1_1_red<-reduce_intered2(intered_X = trs350_1_1_inted, ratio = 0.99)
trs350_1_1_interd<-rbind(torus350_colle_set[[1]][[1]][["noizyX"]], trs350_1_1_red)
