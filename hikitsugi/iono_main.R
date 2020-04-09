library(pracma)

bootstrap.evaluate <- function(X,maxdim,maxscale){
  
}

timeseq <- seq(10,to = 100,by = 10)

bstrp.times.test <- lapply(timeseq, function(t){bootstrap.homology(iono,2,2,times = t)})
bstrp.times <- sapply(bstrp.times.test,function(bl){as.numeric(bl[2:5])})
dimnames(bstrp.times) <- list(c("dim1mhole","dim1dhole","dim2mhole","dim2dhole"),timeseq)

bstrp.times.200s.test <- lapply(seq(10,to = 100,by = 10), function(t){bootstrap.homology(iono,2,2,sample.size = 200,times = t)})
bstrp.times.200s <- sapply(bstrp.times.200s.test,function(bl){as.numeric(bl[2:5])})
dimnames(bstrp.times.200s) <- list(c("dim1mhole","dim1dhole","dim2mhole","dim2dhole"),timeseq)

bstrp.times.300s.test <- lapply(seq(10,to = 100,by = 10), function(t){bootstrap.homology(iono,2,2,sample.size = 300,times = t)})
bstrp.times.300s <- sapply(bstrp.times.300s.test,function(bl){as.numeric(bl[2:5])})
dimnames(bstrp.times.300s) <- list(c("dim1mhole","dim1dhole","dim2mhole","dim2dhole"),timeseq)

bstrp.times.349s.test <- lapply(seq(10,to = 100,by = 10), function(t){bootstrap.homology(iono,2,2,sample.size = 349,times = t)})
bstrp.times.349s <- sapply(bstrp.times.349s.test,function(bl){as.numeric(bl[2:5])})
dimnames(bstrp.times.349s) <- list(c("dim1mhole","dim1dhole","dim2mhole","dim2dhole"),timeseq)


iono.subsample.hd <- sapply(1:nrow(iono.mat),function(i)hausdorff_dist(iono.mat,iono.mat[-i,]))

iono.hd.density <- density(iono.subsample.hd)
sum(diff(iono.hd.density$x),iono.hd.density$y[1:(length(iono.hd.density$y)-1)])
sum(diff(iono.hd.density$x),iono.hd.density$y[2:length(iono.hd.density$y)])

show.phom.with.band <- function(band,diag=diagram){
  diag.mat <- diag[[1]]
  max.scale <- attr(diag.mat,"scale")[2]
  polygon(c(max.scale+1,max.scale-h+1,-1,-1),c(max.scale+1,max.scale+1,h-1,0-1),col = "pink",border = 0)
  par(new=T)
  ShowPHOM(diag = diag)
}



l_b <- function(t,haus.dist){
  sum(length(which(haus.dist>t)))/length(haus.dist)
}