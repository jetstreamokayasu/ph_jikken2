library(parallel)


sphere.collect12<-lapply(1:100, function(i){
  sphere <- sphereUnif(500, 2, 1)
  #cat("nsample=", nsample, "data", i, "\n")
  return(list(nsample = 500, noizyX = sphere, diag = 0))
})
save(sphere.collect12, file = "./data/sphere.collect12")
plot3d(sphere.collect12[[1]][["noizyX"]])

sphere.aggr12<-proposedMethodOnly(X=sphere.collect12, maxdim = 2, maxscale = 3, samples = 10)

#parallelテスト
circles <- lapply(1:12,function(x)circleUnif(n = 200))
cl <- makeCluster(2)
parl<-parLapply(cl,circles,ripsDiag,1,1)
stopCluster(cl)
noparl<-lapply(circles,ripsDiag,1,1)


sphere.aggr12.test<-proposedMethodOnly.parallel(X=sphere.collect12[1], maxdim = 2, maxscale = 3, samples = 2)
sphere.aggr12<-proposedMethodOnly.parallel(X=sphere.collect12, maxdim = 2, maxscale = 3, samples = 10)


ahiru.collect<-list(list(nsample = nrow(ahiru.mat), noizyX = ahiru.mat, diag = 0))
ahiru.aggr<-homMethodsComp2compari3(ahiru.collect, 1, 10, 10)

ahiru.dim1<-cyclenumber(ahiru.aggr[[1]])


#ノイズトーラスKDEによるサイクル推定
collect.aggr<-homMethodsKDE(X = data.collect, maxdim = 2)
collect.dim1<-cyclenumber(collect.aggr[[1]], compare=F)
collect.dim2<-cyclenumber(collect.aggr[[2]], compare=F)

collect.arrange<-plotAggr(collect.aggr, data.collect, 2, 1, compare=F, capture="KDE")

collect.dom1<-plotAggr(aggr = hole.aggr2, collect = data.collect, correct = 2, dim = 1)
collect.dim2<-plotAggr(aggr = hole.aggr2, collect = data.collect, correct = 1, dim = 2)


########アヒルが環状データとはどういうことか
ahiru.pca<-prcomp(ahiru.mat)
plot3d(ahiru.pca[["x"]][,1:3])
text3d(ahiru.pca[["x"]][,1:3], texts = 1:72)
plot(ahiru.pca[["x"]][,1:2])
rgl.snapshot("./data/ahiru3d_2.png")

plot(ahiru.pca[["x"]][,1:2], col=rainbow(12)[(row(ahiru.pca[["x"]][, 1:2])[, 1]-1) %/% 6 + 1], cex=1.5, pch=16)

# for(i in 0:11){
#   
#   points(ahiru.pca[["x"]][(i*6):((i+1)*6), 1:2], col=rainbow(12)[i+1], cex=1.5, pch=16)
#   debugText(i, rainbow(12)[i])
# }
rgl.snapshot("./data/ahiru3d_3.png")

plot3d(ahiru.pca[["x"]][,1:3], type="n")
spheres3d(ahiru.pca[["x"]][, 1:3], col=rainbow(12)[(row(ahiru.pca[["x"]][, 1:3])[, 1]-1) %/% 6 + 1] , radius = 0.3)

for(i in 0:11){
  
  #points3d(ahiru.pca[["x"]][(i*6):((i+1)*6), 1:3], col=rainbow(12)[i+1], size=10, pch=16)
  spheres3d(ahiru.pca[["x"]][(i*6):((i+1)*6), 1:3], col=rainbow(12)[i+1], radius = 0.3)

}

#PCAの逆写像の際に平均点ではなく中心点を足した場合の図示
rad<-seq(0, 2*pi, length=20)
plot(cos(rad[6:14]), sin(rad[6:14]), pch=16, cex=2, xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5))
points(mean(cos(rad[6:14])), mean(sin(rad[6:14])), col=2, pch=16, cex=2)
points(cos(rad[10]), sin(rad[10]), col=4, pch=16, cex=2)

tst_pca=prcomp(cbind(cos(rad[6:14]), sin(rad[6:14])))

#平均点を使った場合の逆写像
points(tst_pca[["x"]][,1]%*%t(tst_pca[["rotation"]][,1])+matrix(rep(tst_pca[["center"]], times=9), 9, 2,byrow=T), col=3, pch=16)
abline((-tst_pca[["rotation"]][2,1]/tst_pca[["rotation"]][1,1])*tst_pca[["center"]][1]+tst_pca[["center"]][2], tst_pca[["rotation"]][2,1]/tst_pca[["rotation"]][1,1], lwd=2, col=3)

#中心点を使った場合の逆写像
points(tst_pca[["x"]][,1]%*%t(tst_pca[["rotation"]][,1])+matrix(rep(c(cos(rad[10]), sin(rad[10])), times=9), 9, 2,byrow=T), col=3, pch=16)
abline((-tst_pca[["rotation"]][2,1]/tst_pca[["rotation"]][1,1])*cos(rad[10])+sin(rad[10]), tst_pca[["rotation"]][2,1]/tst_pca[["rotation"]][1,1], lwd=2, col=3)


theta<-seq(0, 2*pi, length=72)



set.seed(102)
nois<-matrix(rnorm(72*2, 0, 0.1), nrow = 72, ncol = 2)


library(RDRToolbox)
d_isomap <- Isomap(data=d_sr, dims=2, k=5)
plotDR(data=d_isomap[[1]], labels=labels, col=color)


library(vegan)
ahiru_dist<-dist(ahiru.mat)

ahiru_iso <- vegan::isomap(ahiru_dist, ndim = 2, k=5)
plot(ahiru_iso[[1]], cex=1.5, pch=16, col=rainbow(12)[row(ahiru_iso[[1]])[, 1] %/% 6 + 1 ])
