require(phacm)
require(tidyverse)
require(TDA)
require(myfs)
require(rgl)
library(seephacm)

#補間前後のPD比較
#300点トーラス

trs15_300_1_1_pd<-ripsDiag(torus15.300subs[[1]][[1]][["noizyX"]], 2, 3, printProgress = T)
plot(trs15_300_1_1_pd[[1]], main="data1")

trs15_in300_1_1_pd<-ripsDiag(torus15.300insubs[[1]][[1]][["noizyX"]], 2, 3, printProgress = T)
plot(trs15_in300_1_1_pd[[1]])

trs15_300_1_2_10_pd<-lapply(2:10, function(k)ripsDiag(torus15.300subs[[1]][[k]][["noizyX"]], 2, 3, printProgress = T))
plot(trs15_300_1_2_10_pd[[1]][[1]])
save2Rdata(trs15_300_1_2_10_pd)

for (k in 1:9) {
  plot(trs15_300_1_2_10_pd[[k]][[1]], main=paste0("data", k+1))
}

trs15_in300_1_1_10_pd<-lapply(1:10, function(k)ripsDiag(torus15.300insubs[[1]][[k]][["noizyX"]], 2, 3, printProgress = T))
save2Rdata(trs15_in300_1_1_10_pd)
plot(trs15_in300_1_1_10_pd[[2]][[1]])

for (k in 1:10) {
  plot(trs15_in300_1_1_10_pd[[k]][[1]], main=paste0("in_data", k))
}

plotPDs(trs15_300_1_2_10_pd)
plotPDs(trs15_in300_1_1_10_pd)

trs15_300_1_2_pl<-calcLandscape(trs15_300_1_2_10_pd[[1]])
trs15_in300_1_2_pl<-calcLandscape(trs15_in300_1_1_10_pd[[2]])


#補間前に不正解なのと補間後に正解したデータセット比較
wrong1_10<-which(torus15.300subs.aggrs[[1]][[2]] < 0.5)[1:10]
trs15_300_1_w1_10_pd<-lapply(wrong1_10, function(k)ripsDiag(torus15.300subs[[1]][[k]][["noizyX"]], 2, 3, printProgress = T))
save2Rdata(trs15_300_1_w1_10_pd)
plotPDs(trs15_300_1_w1_10_pd)

trs15_in300_1_w1_10_pd<-lapply(wrong1_10, function(k)ripsDiag(torus15.300insubs[[1]][[k]][["noizyX"]], 2, 3, printProgress = T))
save2Rdata(trs15_in300_1_w1_10_pd)
plotPDs(trs15_in300_1_w1_10_pd)

trs15_300_1_w1_pl<-calcLandscape(trs15_300_1_w1_10_pd[[1]])
trs15_in300_1_w1_pl<-calcLandscape(trs15_in300_1_w1_10_pd[[1]])

trs15_300_1_w4_pl<-calcLandscape(trs15_300_1_w1_10_pd[[4]])
trs15_in300_1_w4_pl<-calcLandscape(trs15_in300_1_w1_10_pd[[4]])

trs15_300_1_w9_pl<-calcLandscape(trs15_300_1_w1_10_pd[[9]])
trs15_in300_1_w9_pl<-calcLandscape(trs15_in300_1_w1_10_pd[[9]])

trs15_300_1_w1_10_pls<-lapply(1:10, function(k)calcLandscape(trs15_300_1_w1_10_pd[[k]]))
plot_2ndpls(trs15_300_1_w1_10_pls, vert = F)

trs15_in300_1_w1_10_pls<-lapply(1:10, function(k)calcLandscape(trs15_in300_1_w1_10_pd[[k]]))
plot_2ndpls(trs15_in300_1_w1_10_pls)

#サブサンプルのPDを確かめる
torus15.300subs1_2_subs<-lapply(1:10, function(k){
  data<-torus15.300subs[[1]][[2]][["noizyX"]][sample(300, 300*0.8),]
  return(data)
})

torus15.300subs1_2_subs_pds<-lapply(1:10, function(k)
  ripsDiag(torus15.300subs1_2_subs[[k]], maxdimension = 2, maxscale = 3, printProgress = T))

plotPDs(torus15.300subs1_2_subs_pds)
save2Rdata(torus15.300subs1_2_subs_pds)

torus15.300subs1_2_subs_pls1_5<-plot_pls(torus15.300subs1_2_subs_pds[1:5])
torus15.300subs1_2_subs_pls6_10<-plot_pls(torus15.300subs1_2_subs_pds[6:10])

torus15.300subs1_2_subs_pls<-lapply(1:10, function(k)calcLandscape(torus15.300subs1_2_subs_pds[[k]]))

plot_2ndpls(torus15.300subs1_2_subs_pls)

torus15.300insubs1_2_subs<-lapply(1:10, function(k){
  npoints<-nrow(torus15.300insubs[[1]][[2]][["noizyX"]])
  data<-torus15.300insubs[[1]][[2]][["noizyX"]][sample(npoints, npoints*0.8),]
  return(data)
})

torus15.300insubs1_2_subs_pds<-lapply(1:10, function(k)
  ripsDiag(torus15.300insubs1_2_subs[[k]], maxdimension = 2, maxscale = 3, printProgress = T))

save2Rdata(torus15.300insubs1_2_subs_pds)
plotPDs(torus15.300insubs1_2_subs_pds)

torus15.300insubs1_2_subs_pls<-lapply(1:10, function(k)calcLandscape(torus15.300insubs1_2_subs_pds[[k]]))
plot_2ndpls(torus15.300insubs1_2_subs_pls)

#うまくいかない理由を調べる

inwrong<-which(torus15.300insubs1.aggr[[2]] < 0.5 | torus15.300insubs1.aggr[[2]]>=1.5)

trs15_300_1_iw_pd<-lapply(inwrong, function(k)ripsDiag(torus15.300subs[[1]][[k]][["noizyX"]], 2, 3, printProgress = T))
save2Rdata(trs15_300_1_iw_pd)

plotPDs(trs15_300_1_iw_pd)

trs15_300_1_iw_pls<-lapply(1:3, function(k)calcLandscape(trs15_300_1_iw_pd[[k]]))
plot_2ndpls(trs15_300_1_iw_pls)

trs15_in300_1_iw_pd<-lapply(inwrong, function(k)ripsDiag(torus15.300insubs[[1]][[k]][["noizyX"]], 2, 3, printProgress = T))
save2Rdata(trs15_in300_1_iw_pd)

plotPDs(trs15_in300_1_iw_pd)

trs15_in300_1_iw_pls<-lapply(1:3, function(k)calcLandscape(trs15_in300_1_iw_pd[[k]]))
plot_2ndpls(trs15_in300_1_iw_pls)

#補間点の誤差を調べる
trs15_in300_1_er<-calc_error(torus15.300insubs1_3[[1]][[1]][["noizyX"]], maxr = 2.5, minr = 1, nps = 300)

trs15_in300_1ers<-lapply(1:100, function(i)calc_error(torus15.300insubs1_3[[1]][[i]][["noizyX"]], maxr = 2.5, minr = 1, nps = 300))
par(mgp=c(2.5,1,0))
boxplot(trs15_in300_1ers[1:50], xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6)

oldpar <- par(no.readonly=T)

trs15_in300_1_1_der<-torus_disterror(torus15.300insubs1_3[[1]][[1]][["noizyX"]], maxr = 2.5, minr = 1, nps = 300)
hist(trs15_in300_1_1_der, col="#993435")


trs15_in300_1_w1_10ers<-lapply(wrong1_10, function(i)torus_disterror(torus15.300insubs1_3[[1]][[i]][["noizyX"]], maxr = 2.5, minr = 1, nps = 300))
par(mgp=c(2.5,1,0))
boxplot(trs15_in300_1_w1_10ers, xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6)

#PL一括表示
oldpar<-par(no.readonly=T)
par(cex.lab=3, cex.main=3, cex.axis=3, plt = c(0.2, 0.9, 0.2, 0.9), mfrow=c(2, 5), mgp=c(3.5,1,0))

bpls<-lapply(trs15_300_1_w1_10_pls, function(land){
  
  plot(land[[1]], land[["Land.dim2"]], type = "l", col=3, xlab = "(Birth + Death) / 2", ylab = "(Death - Birth) / 2", ylim=c(0, 0.3))
  abline(h=land[["thresh"]]/2)
  
})
save2Rdata(trs15_300_1_w1_10_pls)

#mtext(side = 3, line=1, outer=T, text = "Title", cex=2)

apls<-lapply(trs15_in300_1_w1_10_pls, function(land){
  
  plot(land[[1]], land[["Land.dim2"]], type = "l", col=3, xlab = "(Birth + Death) / 2", ylab = "(Death - Birth) / 2", ylim=c(0, 0.3))
  abline(h=land[["thresh"]]/2)
  
})
save2Rdata(trs15_in300_1_w1_10_pls)

bspls<-lapply(torus15.300subs1_2_subs_pls, function(land){
  
  plot(land[[1]], land[["Land.dim2"]], type = "l", col=3, xlab = "(Birth + Death) / 2", ylab = "(Death - Birth) / 2", ylim=c(0, 0.3))
  abline(h=land[["thresh"]]/2)
  
})
save2Rdata(torus15.300subs1_2_subs_pls)

aspls<-lapply(torus15.300insubs1_2_subs_pls, function(land){
  
  plot(land[[1]], land[["Land.dim2"]], type = "l", col=3, xlab = "(Birth + Death) / 2", ylab = "(Death - Birth) / 2", ylim=c(0, 0.3))
  abline(h=land[["thresh"]]/2)
  
})
save2Rdata(torus15.300insubs1_2_subs_pls)

par(oldpar)

#補間前後のトーラスプロット
figurePlot(torus15.300subs[[1]][[2]][["noizyX"]][1:300,])
rgl.postscript("./data/b_torus.eps", fmt="eps" ) 

points3d(torus15.300insubs[[1]][[2]][["noizyX"]][300:torus15.300insubs[[1]][[2]][["nsample"]] ,], col=2)
rgl.postscript("./data/a_torus.eps", fmt="eps" ) 

trs15_in300_1ders<-lapply(1:100, function(i)torus_disterror(torus15.300insubs1_3[[1]][[i]][["noizyX"]], maxr = 2.5, minr = 1, nps = 300))
boxplot(trs15_in300_1ders[1:50], xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6)

#animation試し
library(animation)
library(jpeg)
library(magick)
library("tcltk2")

#基本的な使い方:saveGIFコマンド
#画像の切替間隔(秒):intervalオプション
#ファイル名:movie.nameオプション
#作業フォルダに出力されます
saveGIF({
  #繰り返し回数を設定する 
  for (i in 1:5){
    #処理内容を記述する
    plot(runif(10), ylim = c(0,1))
  }
}, interval = 1.0, movie.name = "TEST.gif")

###「jpeg」パッケージと連携して写真をGIFアニメーション化#####
#jpgファイルの保存フォルダを指定
OpenFiles <- paste(as.character(tkchooseDirectory(title = "~/R/実験データ/2019_02_10"), sep = "", collapse =""))

#OpenFiles<-todaywd()

#保存フォルダを作業フォルダに変更
setwd(OpenFiles)

#フォルダ内のファイルを取得
ItemList <- list.files(path = OpenFiles)

#gifアニメーションを作成
par(mar = c(0, 0, 0, 0))
par(oma = c(0, 0, 0, 0))
saveGIF(
  
  for (i in seq(length(ItemList))){
    
    PlotJpg <- readJPEG(ItemList[i])
    plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
    rasterImage(PlotJpg, 0, 1, 1, 0)
    
  },
  
  interval = 1.0, movie.name = "Animetion.gif")

wave <- function() {
  for(t in 1:100) {
    plot(function(x){ sin(x + 0.08 * pi * t) },
         -pi, 2*pi, xlab="x", ylab="sin(x)",
         col="blue", lwd=3)
  }
}

saveGIF(wave(), interval=0.05, moviename="wave",
          movietype="gif", outdir=getwd(),
          width=640, height=480)

xa<-c(-3, -2, -1, 0)
yb<-c(2, 0, 3, 1)
plot(xa, yb, pch=20, xlim=c(-6, 4), ylim=c(-3, 7))

plot.circle <- function(x, y, r){
  theta <- seq(-pi, pi, length=100)
  polygon(x + r*cos(theta), y + r*sin(theta), col="gray")
}

sapply(1:6, function(k)plot.circle(xa[k], yb[k], sqrt(17)/2))

saveGIF({
  rs<-seq(0.1, 2, length=20)
  cat("alpha\n")
  for (i in 1:20){
    r<-rs[i]
    theta <- seq(-pi, pi, length=100)
    cat("r=", r)
    plot(xa, yb, pch=20, xlim=c(-6, 4), ylim=c(-3, 7), typ="n")
    for(k in 1:4){
      polygon(xa[k] + r*cos(theta), yb[k] + r*sin(theta), col="gray")
      #debugText(xa[k], yb[k])
      }
    #plot.circle(xa[i], yb[i], 0.5)
    par(new=T)
    plot(xa, yb, pch=20, xlim=c(-6, 4), ylim=c(-3, 7))
  }
}, interval = 0.5, movie.name = "circ.gif", width=640, height=640)



xa<-c(-3, -2, -1, 0)
yb<-c(2, 0, 3, 0.5)
xy<-cbind(xa, yb)
xy.dist<-dist(xy) %>% as.matrix()
xydist<-dist(xy)

oopts <- ani.options(ffmpeg = paste0(getwd(), "//ffmpeg.exe"), interval = 0.1, nmax=1000)
ani.options(oopts)
saveVideo({
  rs<-seq(0.05, 3.4/2, length=100)
  theta <- seq(-pi, pi, length=100)
  cat("start\n")
  for (i in 1:100){
    r<-rs[i]
    plot(xa, yb, xlim=c(-6, 4), ylim=c(-3, 7), typ="n")
    for(k in 1:4){
      polygon(xa[k] + r*cos(theta), yb[k] + r*sin(theta), col="pink")
      
    }
    #cat("TF=", (xy.dist < r & xy.dist > 0), "\n")
    #cat("r=", r, "\n")
    if(T %in% (xy.dist < 2*r & xy.dist > 0)){
      
      grp<-which(xy.dist!=0 & xy.dist < r*2, arr.ind=TRUE)
      #cat("dist=", grp, "\n")
      # cat("nrow=", nrow(grp), "\n")
      #  cat("xy1=", xy[grp[j, 1], ], "\n")
      #  cat("xy2=", xy[grp[j, 2], ], "\n")
      # cat("r=", r, "\n")
      for(j in 1:nrow(grp))
      lines(c(xy[grp[j, 1], 1], xy[grp[j, 2], 1]), c(xy[grp[j, 1], 2], xy[grp[j, 2],2]), col=2)
      # cat("nrow=", nrow(grp), "\n")
      # cat("xy1=", xy[grp[j, 1], ], "\n")
      # cat("xy2=", xy[grp[j, 2], ], "\n")
      #cat("r=", r, "\n")
    }
    
    par(new=T)
    plot(xa, yb, pch=20, xlim=c(-6, 4), ylim=c(-3, 7), cex=2)
    
    if(r*2 > xydist[order(-xydist)[2]]){
      
      polygon(c(-2, 0, -1, -3), c(0, 0.5, 3, 2), density = 10, angle = -45, col=1)
      
    }
  }
}, video.name = "circ.mp4")


saveGIF({
  rs<-seq(0.05, 3.4/2, length=100)
  theta <- seq(-pi, pi, length=100)
  cat("start\n")
  for (i in 1:100){
    r<-rs[i]
    plot(xa, yb, xlim=c(-6, 4), ylim=c(-3, 7), typ="n")
    for(k in 1:4){
      polygon(xa[k] + r*cos(theta), yb[k] + r*sin(theta), col="pink")
      
    }
    #cat("TF=", (xy.dist < r & xy.dist > 0), "\n")
    #cat("r=", r, "\n")
    if(T %in% (xy.dist < 2*r & xy.dist > 0)){
      
      grp<-which(xy.dist!=0 & xy.dist < r*2, arr.ind=TRUE)
      #cat("dist=", grp, "\n")
      # cat("nrow=", nrow(grp), "\n")
      #  cat("xy1=", xy[grp[j, 1], ], "\n")
      #  cat("xy2=", xy[grp[j, 2], ], "\n")
      # cat("r=", r, "\n")
      for(j in 1:nrow(grp))
        lines(c(xy[grp[j, 1], 1], xy[grp[j, 2], 1]), c(xy[grp[j, 1], 2], xy[grp[j, 2],2]), col=2)
      # cat("nrow=", nrow(grp), "\n")
      # cat("xy1=", xy[grp[j, 1], ], "\n")
      # cat("xy2=", xy[grp[j, 2], ], "\n")
      #cat("r=", r, "\n")
    }
    
    par(new=T)
    plot(xa, yb, pch=20, xlim=c(-6, 4), ylim=c(-3, 7), cex=2)
    
    if(r*2 > xydist[order(-xydist)[2]]){
      
      polygon(c(-2, 0, -1, -3), c(0, 0.5, 3, 2), density = 10, angle = -45, col=1)
      
    }
  }
}, interval = 0.1, movie.name = "test.gif")



#フィルトレーション画像の作成
xa<-c(-3, -2, -1, 0)
yb<-c(2, 0, 3, 0.5)
xy<-cbind(xa, yb)
xy_dist<-dist(xy) %>% as.matrix()
xydist<-dist(xy)
#1.点群
plot(xa, yb, pch=20, xlim=c(-6, 4), ylim=c(-3, 7), cex=2, bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

#2.フィルトレーション開始
for(k in 1:4){
     polygon(xa[k] + cos(theta), yb[k] + sin(theta), col="pink")
      
}
par(new=T)
plot(xa, yb, pch=20, xlim=c(-6, 4), ylim=c(-3, 7), cex=2, bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

#3.穴の発生
r<-xydist[order(-xydist)[3]]/2
for(k in 1:4){
  polygon(xa[k] + r*cos(theta), yb[k] + r*sin(theta), col="pink")
  
}

par(new=T)
plot(xa, yb, pch=20, xlim=c(-6, 4), ylim=c(-3, 7), cex=2, bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

grp<-which(xy.dist!=0 & xy.dist < r*2, arr.ind=TRUE)
for(j in 1:nrow(grp)){
  lines(c(xy[grp[j, 1], 1], xy[grp[j, 2], 1]), c(xy[grp[j, 1], 2], xy[grp[j, 2],2]))}
lines(c(0, -1), c(0.5, 3))

#4.穴の消滅
r<-xydist[order(-xydist)[2]]/2
for(k in 1:4){
  polygon(xa[k] + r*cos(theta), yb[k] + r*sin(theta), col="pink")
  
}

par(new=T)
plot(xa, yb, pch=20, xlim=c(-6, 4), ylim=c(-3, 7), cex=2, bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

grp<-which(xy.dist!=0 & xy.dist < r*2, arr.ind=TRUE)
for(j in 1:nrow(grp)){
  lines(c(xy[grp[j, 1], 1], xy[grp[j, 2], 1]), c(xy[grp[j, 1], 2], xy[grp[j, 2],2]))}
lines(c(-2, -1), c(0, 3))
polygon(c(-2, 0, -1, -3), c(0, 0.5, 3, 2), density = 10, angle = -45, col=1)
