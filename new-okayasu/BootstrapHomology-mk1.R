bootstrap.homology <- function(X,maxdim,maxscale,const.band=0,maximum.thresh = F){
  require(TDA)
  # require(pracma)
  if(!("bootsSamples" %in% class(X))) stop("input must be bootsSamples")
  peak <- matrix(0,maxdim,length(X))
  # band <- ifelse(const.band > 0,const.band,hausdInterval(X, m=sample.size, B=times, alpha = (1-confidence)))
  tseq <- seq(0,maxscale,length.out = 1000)
  diags <- lapply(X,function(x)calcPhom(x,maxdim,maxscale,ret = T,plot = F))
  print(sapply(diags,function(diag)calcDiagCentroid(diag)[3]))
  band <- ifelse(const.band==0,max(sapply(diags,function(diag)calcDiagCentroid(diag)[3])),const.band)
  print(band)
  
  for (t in 1:length(X)) {
    land <- lapply(1:maxdim,function(d)landscape(diags[[t]][[1]],dimension = d,KK = 1,tseq = tseq))
    if(maximum.thresh) band <- max(sapply(land,max))/4
    for(d in 1:maxdim){
      peak[d,t] <- calc.landscape.peak(X=land[[d]], thresh = (band/(2*d)), tseq=tseq)
    }
  }
  
  dimnames(peak) <- list(paste0("dim",1:maxdim),paste0("sample",1:length(X)))
  bootstrap.summary <- list(peak=peak)
  bootstrap.summary <- append(bootstrap.summary,c(band=band,show.hole.density(peak)))
  class(bootstrap.summary) <- "smoothPhom"
  return(bootstrap.summary)
}

calcDiagCentroid <- function(diag = diagram){
  if(class(diag)=="list") diag <- diag[[1]]
  diag <- diag[-which(diag[,1]==0),]
  centroid <- apply(diag[,2:3],2,mean)
  cpersistence <- (centroid[2]-centroid[1])
  ret <- c(centroid,cpersistence,cpersistence*2)
  names(ret) <- c("birth","death","persistence","noizes.thresh")
  return(ret)
}

plot.smoothPhom <- function(x) x <- show.hole.density(x$peak)

#平滑化して平均ピーク数を数える。ダイアグラム直入力可。
calc.landscape.peak <- function(X = diagram,dimension = 1,spar = seq(0,1,0.1),thresh = 0,tseq = tseq,show = F){
  if (class(X)=="list") {
    require(TDA)
    X <- X[["diagram"]]
    tseq  <- seq(0,attr(X,"scale")[2],length.out = 1000)
    X <- landscape(X,dimension = dimension,1,tseq)
  }
  else if(missing(tseq)) tseq <- 1:length(X)
  
  smoothed.land <- lapply(spar, function(sp){smooth.spline(x=tseq,y=X,spar=sp)})
  if(show){ lapply(1:length(spar),function(i) {
    if(i==1) plot(tseq,smoothed.land[[i]]$y,type = "l",col=rainbow(length(spar))[i],ylim = c(0,max(X)),
                  xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2")
    else plot(tseq,smoothed.land[[i]]$y,type = "l",col=rainbow(length(spar))[i],ann=F,axes = F,ylim = c(0,max(X)))
    par(new=T)})
    abline(thresh,0,col=2)}
  # for single showing
  # if(show) lapply(1:length(spar),function(i) plot(tseq,smoothed.land[[i]]$y,type = "l",col=1,ann = F,ylim = c(0,max(X))))
  par(new=F)
  return(mean(sapply(smoothed.land, count.local.maximal, T, thresh)))
}


#1次関数のピーク数を数える Internal Function
count.local.maximal <- function(x,weakcut = T,thresh = max(x)/4,show.thresh=F){
  if (class(x)=="smooth.spline") x <- x$y
  if (thresh==0) thresh <- max(x)/4
  peak <- 0
  d <- diff(x)
  thresh <- thresh * weakcut
  for (i in 1:(length(x)-2)) {
    if (d[i] >= 0 && d[i+1] < 0 && x[i] > thresh) {
      peak <- peak + 1
    }
  }
  if(show.thresh)print(paste("thresh",thresh))
  return(peak)
}

#各次元毎のpeak数行列を密度推定し，表示する。 Internal FUnction
show.hole.density <- function(X){
  dens <- apply(X, 1, density)
  maxdim <- nrow(X)
  bootstrap.summary <- list()
  xlim = c(min(sapply(dens, function(den)min(den$x))),max(sapply(dens, function(den)max(den$x))))
  ylim = c(0,max(sapply(dens, function(den)max(den$y))))
  plot(1,type = "n",xlim = xlim,ylim = ylim,xlab = "betti number",ylab = "probability")
  for (d in 1:maxdim) {
    mhole <- mean(X[d,])
    dhole <- dens[[d]][["x"]][which.max(dens[[d]][["y"]])]
    par(new=T)
    plot(dens[[d]],xlim=xlim,ylim=ylim,col=d+1,ann = F,axes = F)
    print(paste0("dimension ",d,", ",round(mhole,digits = 2)," mean hole, ",round(dhole)," density hole"))
    bootstrap.summary[[paste0("dim",d,"mhole")]] <- mhole
    bootstrap.summary[[paste0("dim",d,"dhole")]] <- dhole
  }
  legend("topright",legend = paste("dim",1:maxdim),pch = 0,col = 2:(maxdim+1))
  return(bootstrap.summary)
}

calc.bootstrap.confidence <- function(confidence,maxscale,haus.dist,digit=3){
  Lb <- function(t,haus.dist) sum(length(which(haus.dist>t)))/length(haus.dist)
  max.dist <- max(haus.dist)
  from <- 0
  to <- max.dist
  for (i in 1:digit) {
    tseq <- seq(from,to = to,length.out = 10)
    # ErrorText(from,to,haus.dist,seq(from,to = to,length.out = 10))
    prob <- sapply(tseq, Lb, haus.dist)
    # print(prob)
    from <- tseq[which(prob < confidence)[1]-1]
    # print(from)
    to <- tseq[which(prob < confidence)[1]]
  }
  c <- 2*from
  return(c)
}

#各次元ごとに閾値を求めている
#ワイブル分布による閾値設定を試す
bootstrap.homology.mk1 <- function(X,maxdim,maxscale, const.band=0, maximum.thresh = F){
  require(TDA)
  # require(pracma)
  if(!("bootsSamples" %in% class(X))) stop("input must be bootsSamples")
  peak <- matrix(0,maxdim,length(X))
  # band <- ifelse(const.band > 0,const.band,hausdInterval(X, m=sample.size, B=times, alpha = (1-confidence)))
  tseq <- seq(0,maxscale,length.out = 1000)
  diags <- lapply(X,function(x)calcPhom(x,maxdim,maxscale,ret = T,plot = F))
  print(sapply(diags,function(diag)threshWeibull(diag)))
  band <- ifelse(const.band==0,max(sapply(diags,function(diag)threshWeibull(diag))),const.band)
  print(band)
  
  for (t in 1:length(X)) {
    land <- lapply(1:maxdim,function(d)landscape(diags[[t]][[1]],dimension = d,KK = 1,tseq = tseq))
    if(maximum.thresh) band <- max(sapply(land,max))/4
    for(d in 1:maxdim){
      peak[d,t] <- calc.landscape.peak(X=land[[d]], thresh = band/(2*d), tseq=tseq)
    }
  }
  
  dimnames(peak) <- list(paste0("dim",1:maxdim),paste0("sample",1:length(X)))
  bootstrap.summary <- list(peak=peak)
  bootstrap.summary <- append(bootstrap.summary,c(band=band,show.hole.density(peak)))
  class(bootstrap.summary) <- "smoothPhom"
  return(bootstrap.summary)
}

#各次元ごとに、その次元の平均を閾値としている
calcDiagCentroid.mk1 <- function(diag = diagram){
  if(class(diag)=="list") diag <- diag[[1]]
  diag <- diag[-which(diag[,1]==0),]
  centroid<-matrix(0, 2, 2)
  centroid[1,] <- apply(diag[diag[,1]==1,2:3],2,mean)
  centroid[2,] <- apply(diag[diag[,1]==2,2:3],2,mean)
  cpersistence <- c(centroid[1,2]-centroid[1,1], centroid[2,2]-centroid[2,1])
  ret <- list(centroid=centroid,cpersistence=cpersistence,noizes.thresh=cpersistence*2)
  #names(ret) <- c("birth","death","persistence","noizes.thresh")
  return(ret)
}

#2次のパーシステンスを2倍にしてから平均を求める
calcDiagCentroid.mk2 <- function(diag = diagram){
  if(class(diag)=="list") diag <- diag[[1]]
  diag <- diag[-which(diag[,1]==0),]
  centroid1 <- diag[diag[,1]==1,3]-diag[diag[,1]==1,2]
    centroid2 <- (diag[diag[,1]==2,3]-diag[diag[,1]==2,2])*2
    cpersistence <- mean(c(centroid1, centroid2))
  ret <- c(cpersistence,noizes.thresh=cpersistence*2)
  names(ret) <- c("cpersistence","noizes_thresh")
  return(ret)
}

#2次のパーシステンスを2倍にしてから求めた平均を閾値としている
bootstrap.homology.mk2 <- function(X,maxdim,maxscale,const.band=0,maximum.thresh = F){
  require(TDA)
  # require(pracma)
  if(!("bootsSamples" %in% class(X))) stop("input must be bootsSamples")
  peak <- matrix(0,maxdim,length(X))
  # band <- ifelse(const.band > 0,const.band,hausdInterval(X, m=sample.size, B=times, alpha = (1-confidence)))
  tseq <- seq(0,maxscale,length.out = 1000)
  diags <- lapply(X,function(x)calcPhom(x,maxdim,maxscale,ret = T,plot = F))
  print(sapply(diags,function(diag)calcDiagCentroid.mk2(diag)[1]))
  band <- ifelse(const.band==0,max(sapply(diags,function(diag)calcDiagCentroid.mk2(diag)[1])),const.band)
  print(band)
  
  for (t in 1:length(X)) {
    land <- lapply(1:maxdim,function(d)landscape(diags[[t]][[1]],dimension = d,KK = 1,tseq = tseq))
    if(maximum.thresh) band <- max(sapply(land,max))/4
    for(d in 1:maxdim){
      peak[d,t] <- calc.landscape.peak(X=land[[d]], thresh = (band/d), tseq=tseq)
    }
  }
  
  dimnames(peak) <- list(paste0("dim",1:maxdim),paste0("sample",1:length(X)))
  bootstrap.summary <- list(peak=peak)
  bootstrap.summary <- append(bootstrap.summary,c(band=band,show.hole.density(peak)))
  class(bootstrap.summary) <- "smoothPhom"
  return(bootstrap.summary)
}

#2次のパーシステンスを2倍にしてから求めた平均を閾値としている
#かつ並行処理
bootstrap.homology.parallel <- function(X,maxdim,maxscale,const.band=0,maximum.thresh = F){
  require(TDA)
  # require(pracma)
  cluster<-makeCluster(2)
  clusterEvalQ(cluster, library(TDA))
  clusterExport(cluster, "calcPhom")
  if(!("bootsSamples" %in% class(X))) stop("input must be bootsSamples")
  peak <- matrix(0,maxdim,length(X))
  # band <- ifelse(const.band > 0,const.band,hausdInterval(X, m=sample.size, B=times, alpha = (1-confidence)))
  tseq <- seq(0,maxscale,length.out = 1000)
  diags <- parLapply(cluster, X, function(x)calcPhom(x, maxdim, maxscale, ret = T, plot = F))
  print(sapply(diags,function(diag)calcDiagCentroid.mk2(diag)[1]))
  band <- ifelse(const.band==0,max(sapply(diags,function(diag)calcDiagCentroid.mk2(diag)[1])),const.band)
  print(band)
  
  for (t in 1:length(X)) {
    land <- lapply(1:maxdim,function(d)landscape(diags[[t]][[1]],dimension = d,KK = 1,tseq = tseq))
    if(maximum.thresh) band <- max(sapply(land,max))/4
    for(d in 1:maxdim){
      
      per.dim0<-diags[[t]][[1]][diags[[t]][[1]][,1]==0, ]
      per.dim0<-per.dim0[-1]
      per.dimd<-calcper(diags[[t]][[1]], d)
      per.dimd<-per.dimd*d
      
      debugText(max(per.dimd), mean(per.dim0))
      
      if(max(per.dimd) < mean(per.dim0)){
        peak[d,t]<-0
      }
     else{ peak[d,t] <- calc.landscape.peak(X=land[[d]], thresh = (band/d), tseq=tseq)}
    }
  }
  
  dimnames(peak) <- list(paste0("dim",1:maxdim),paste0("sample",1:length(X)))
  bootstrap.summary <- list(peak=peak)
  bootstrap.summary <- append(bootstrap.summary,c(band=band,show.hole.density(peak)))
  class(bootstrap.summary) <- "smoothPhom"
  stopCluster(cluster)
  return(bootstrap.summary)
}
