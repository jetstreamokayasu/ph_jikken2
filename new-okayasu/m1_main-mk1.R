require(myfs)

# mean landscape を描画します。
meanLandscape <- function(X,maxdim,maxscale,alpha = 0.95,col=c("red","pink","green","lightgreen"),ret=F){
  if(!("bootsSamples" %in% class(X))) stop("argument X must be bootsSamples")
  diags <- lapply(X,function(x)calcPhom(x,maxdim,maxscale,ret = T))
  lands <- lapply(1:maxdim,function(d,diags)
    lapply(diags,function(diag,d)
      landscape(diag[[1]],dimension = d),d),diags)
  meanland <- lapply(lands,function(dimlands)apply(matlist(dimlands),2,mean))
  bands <- sapply(1:maxdim,function(d)
    meanLandscapeConfidence(lands[[d]], 1, length(X),alpha = alpha))

  obj <- lapply(1:maxdim,function(d){
    attr(meanland[[d]],"band") <- bands[d]
    return(meanland[[d]])
  })
  attr(obj,"col") <- col
  class(obj) <- "meanland"
  
  plot(obj)
  if(ret) return(obj)
}

plot.meanland <- function(x){
  meanland <- x
  col <- attr(meanland,"col")
  plot(meanland[[which.max(sapply(meanland,max))]],type = "n",xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2")
  for(d in 1:length(meanland)){
    band <- attr(meanland[[d]],"band")
    polygon(c(1:500,500:1),c(meanland[[d]]+band,rev(meanland[[d]]-band)),
            border = col[(d-1)*2+1],col = col[(d-1)*2+2])
    lines(meanland[[d]])
  }
}

# mean landscape の信頼区間を計算します。Internal Function
meanLandscapeConfidence <- function(lands, dim, boots.samples, alpha){
  meanland <- apply(matlist(lands),2,mean)
  theta = numeric(boots.samples)
  for(i in 1:boots.samples){
    xi <- rnorm(boots.samples)
    theta[i] = max(length(lands)^(-1/2)*abs(sum(xi*(lands[[i]]-meanland))))
  }
  Z = min(which(sapply(seq(0,1,by = 0.001),function(z)mean(theta > z))<=alpha))*0.001
  return(Z)
}

# persistent diagramのホモロジー類を数えます。
calcDiagHomology <- function(thresh,diag=diagram){
  if(class(diag)=="list") diag <- diag[[1]]
  maxscale <- attr(diag,"scale")[2]
  maxdim <- attr(diag,"maxdimension")
  confHom <- diag[which((diag[,3]-diag[,2])/2 > thresh),]
  confHom <- try(as.matrix(confHom[-which(confHom[,3]==maxscale),]),silent = T)
  if(class(confHom)=="try-error") return(c(0,0,0))
  nhom <- sapply(0:maxdim, function(d) sum(nrow(confHom[which(confHom[,1]==d),])))
  names(nhom) <- paste0("dim",0:maxdim)
  return(nhom)
}

bootstrapper <- function(X,size,samples){
  require(TDA)
  hausdint <- hausdInterval(X = X,m = size,B = samples)
  X <- lapply(1:samples,function(i)X[sample(nrow(X),size),])
  attr(X,"size") <- size
  attr(X,"samples") <- samples
  attr(X,"hausd") <- hausdint
  class(X) <- c("matrix","bootsSamples")#R6クラスに変更できると良い
  return(X)
}

# 4種の手法を比較
#指数分布に基づく閾値を設定
homMethodsComp2compari <- function(X,maxdim,maxscale,samples){
  aggr1 <- matrix(0,length(X),5)
  aggr2 <- matrix(0,length(X),5)
  dimnames(aggr1) <- list(paste0("data-set", 1:length(X)), c("confidence","bestconf","landscape","meanland","proposed"))
  dimnames(aggr2) <- dimnames(aggr1)
  
  for(t in 1:length(X)){
    
    cat("data set", t, "calculating\n")
    
    size<-X[[t]]$nsample*(3/5)
    B <- bootstrapper(X[[t]]$noizyX,size,samples)
    X[[t]]$diag<-calcPhom(X[[t]]$noizyX,maxdim,maxscale, plot = T,ret = T)
    m1 <- calcDiagHomology(attr(B,"hausd"))[2:3]
    m2 <- calcDiagHomology(hausdInterval(X[[t]]$noizyX,nrow(X[[t]]$noizyX)-1))[2:3]
    xland <- lapply(1:maxdim,function(d)landscape(diagram[[1]],d))
    m3 <- sapply(xland,count.local.maximal,T,max(sapply(xland,max))/4)
    mland <- meanLandscape(B,maxdim,maxscale,ret = T)
    m4 <- sapply(1:maxdim,function(d)count.local.maximal(mland[[d]],thresh = max(sapply(xland,max))/4))
    
    band<-threshDetermine(X[[t]]$diag,maxscale)
    speak <- bootstrap.homology.mk2(B,maxdim,maxscale, const.band = band)
    m5 <- sapply(1:maxdim,function(d)speak[[paste0("dim",d,"dhole")]])
    
    aggr1[t,1] <- m1[1]
    aggr1[t,2] <- m2[1]
    aggr1[t,3] <- m3[1]
    aggr1[t,4] <- m4[1]
    aggr1[t,5] <- m5[1]
    
    aggr2[t,1] <- m1[2]
    aggr2[t,2] <- m2[2]
    aggr2[t,3] <- m3[2]
    aggr2[t,4] <- m4[2]
    aggr2[t,5] <- m5[2]
  }
  Xdiag<-lapply(1:length(X), function(k){return(X[[k]]$diag)})
  aggrs <- list(aggr1,aggr2)
  aggrs <- append(aggrs,list(Xsize=nrow(X[[1]]),Xsamples=length(X),
                             Bsize=size,Bsamples=samples,
                             maxdim=maxdim,maxscale=maxscale, Xdiag=Xdiag))
  class(aggrs) <- "bettiComp"
  
  return(aggrs)
}

plot.bettiComp <- function(x, dim=1, vline=2, xlim=c(0,8), ylim=c(1,max(x[[dim]])), emphasis = 1){
  x <- x[[dim]]
  xmax <- xlim[2]
  plot(1,xlim = xlim,ylim = ylim,xlab = "number of holes",ylab = "times")
  abline(v=vline)
  for(i in 1:nrow(x)){
    par(new=T)
    if(i==emphasis){
      plot(0:xmax,x[i,1:(xmax+1)],xlim = xlim,ylim = ylim,col=rainbow(nrow(x))[i],ann = F,axes = F, pch=16)
      lines(0:xmax,x[i,1:(xmax+1)],col=rainbow(nrow(x))[i],lwd=3)
    }
    else if(i != 2){
    plot(0:xmax,x[i,1:(xmax+1)],xlim = xlim,ylim = ylim,col=rainbow(nrow(x))[i],ann = F,axes = F, pch=16)
    lines(0:xmax,x[i,1:(xmax+1)],col=rainbow(nrow(x))[i],lty=i+1)
    }
  }
  legend("topright",legend = rownames(x)[-2],pch = 3, col=rainbow(nrow(x))[-2])
}

homologyMethodsRandomDataComp <- function(Xsamples){
  Xsize <- rnorm(Xsamples,600,100)
  X <- lapply(Xsize,torusUnif,1,2.5)
  maxdim <- 2
  maxscale <- 3
  Bsize <- rnorm(Xsamples,0.6,0.1)
  while(which(Bsize>=1)!=0) Bsize <- rnorm(Xsamples,0.6,0.1)
  Bsize <- round(Bsize*Xsize)
  Bsamples <- 10
    
  aggr1 <- matrix(0,5,10)
  aggr2 <- matrix(0,5,10)
  dimnames(aggr1) <- list(c("confidence","bestconf","landscape","meanland","proposed"),paste0(0:19,"holes"))
  dimnames(aggr2) <- dimnames(aggr1)
  
  for(t in 1:length(X)){
    B <- bootstrapper(X[[t]],Bsize,Bsamples)
    calcPhom(X[[t]],maxdim,maxscale)
    m1 <- calcDiagHomology(attr(B,"hausd"))[2:3]
    m2 <- calcDiagHomology(hausdInterval(X[[t]],nrow(X[[t]])-1))[2:3]
    xland <- lapply(1:maxdim,function(d)landscape(diagram[[1]],d))
    m3 <- sapply(xland,count.local.maximal,T,max(sapply(xland,max))/4)
    mland <- meanLandscape(B,maxdim,maxscale,ret = T)
    m4 <- sapply(1:maxdim,function(d)count.local.maximal(mland[[d]],thresh = max(sapply(xland,max))/4))
    speak <- bootstrap.homology(B,maxdim,maxscale)
    m5 <- sapply(1:maxdim,function(d)round(speak[[paste0("dim",d,"dhole")]]))
    aggr1[1,m1[1]+1] <- aggr1[1,m1[1]+1] + 1
    aggr1[2,m2[1]+1] <- aggr1[2,m2[1]+1] + 1
    aggr1[3,m3[1]+1] <- aggr1[3,m3[1]+1] + 1
    aggr1[4,m4[1]+1] <- aggr1[4,m4[1]+1] + 1
    aggr1[5,m5[1]+1] <- aggr1[5,m5[1]+1] + 1
    
    aggr2[1,m1[2]+1] <- aggr2[1,m1[2]+1] + 1
    aggr2[2,m2[2]+1] <- aggr2[2,m2[2]+1] + 1
    aggr2[3,m3[2]+1] <- aggr2[3,m3[2]+1] + 1
    aggr2[4,m4[2]+1] <- aggr2[4,m4[2]+1] + 1
    aggr2[5,m5[2]+1] <- aggr2[5,m5[2]+1] + 1
  }
  aggrs <- list(aggr1,aggr2)
  aggrs <- append(aggrs,list(Xsize=Xsize,Xsamples=length(X),
                             Bsize=Bsize*Xsize,Bsamples=Bsamples,
                             maxdim=maxdim,maxscale=maxscale))
  class(aggrs) <- "bettiComp"
  
  return(aggrs)
}

#各次元のごとに閾値を求めている
#ワイブル分布による閾値設定を試す
homMethodsComp2compari2 <- function(X,maxdim,maxscale,samples){
  aggr1 <- matrix(0,length(X),5)
  aggr2 <- matrix(0,length(X),5)
  dimnames(aggr1) <- list(paste0("data-set", 1:length(X)), c("confidence","bestconf","landscape","meanland","proposed"))
  dimnames(aggr2) <- dimnames(aggr1)
  
  for(t in 1:length(X)){
    
    cat("data set", t, "calculating\n")
    
    size<-X[[t]]$nsample*(3/5)
    B <- bootstrapper(X[[t]]$noizyX,size,samples)
    X[[t]]$diag<-calcPhom(X[[t]]$noizyX,maxdim,maxscale, plot = T,ret = T)
    m1 <- calcDiagHomology(attr(B,"hausd"))[2:3]
    m2 <- calcDiagHomology(hausdInterval(X[[t]]$noizyX,nrow(X[[t]]$noizyX)-1))[2:3]
    xland <- lapply(1:maxdim,function(d)landscape(diagram[[1]],d))
    m3 <- sapply(xland,count.local.maximal,T,max(sapply(xland,max))/4)
    mland <- meanLandscape(B,maxdim,maxscale,ret = T)
    m4 <- sapply(1:maxdim,function(d)count.local.maximal(mland[[d]],thresh = max(sapply(xland,max))/4))
    
    # ret<-calcDiagCentroid.mk1(X[[t]]$diag)
    # debugText(ret)
    speak <- bootstrap.homology.mk1(B,maxdim,maxscale)
    m5 <- sapply(1:maxdim,function(d)speak[[paste0("dim",d,"dhole")]])
    
    aggr1[t,1] <- m1[1]
    aggr1[t,2] <- m2[1]
    aggr1[t,3] <- m3[1]
    aggr1[t,4] <- m4[1]
    aggr1[t,5] <- m5[1]
    
    aggr2[t,1] <- m1[2]
    aggr2[t,2] <- m2[2]
    aggr2[t,3] <- m3[2]
    aggr2[t,4] <- m4[2]
    aggr2[t,5] <- m5[2]
  }
  Xdiag<-lapply(1:length(X), function(k){return(X[[k]]$diag)})
  aggrs <- list(aggr1,aggr2)
  aggrs <- append(aggrs,list(Xsize=nrow(X[[1]]),Xsamples=length(X),
                             Bsize=size,Bsamples=samples,
                             maxdim=maxdim,maxscale=maxscale, Xdiag=Xdiag))
  class(aggrs) <- "bettiComp"
  
  return(aggrs)
}

#もともとの関数
homMethodsComp2compari3 <- function(X,maxdim,maxscale,samples){
  aggr1 <- matrix(0,length(X),5)
  aggr2 <- matrix(0,length(X),5)
  dimnames(aggr1) <- list(paste0("data-set", 1:length(X)), c("confidence","bestconf","landscape","meanland","proposed"))
  dimnames(aggr2) <- dimnames(aggr1)
  
  for(t in 1:length(X)){
    
    cat("data set", t, "calculating\n")
    
    size<-round(X[[t]]$nsample*(4/5))
    B <- bootstrapper(X[[t]][["noizyX"]],size,samples)
    X[[t]]$diag<-calcPhom(X[[t]]$noizyX,maxdim,maxscale, plot = T,ret = T)
    m1 <- calcDiagHomology(attr(B,"hausd"))[2:3]
    m2 <- calcDiagHomology(hausdInterval(X[[t]]$noizyX,nrow(X[[t]]$noizyX)-1))[2:3]
    xland <- lapply(1:maxdim,function(d)landscape(diagram[[1]],d))
    m3 <- sapply(xland,count.local.maximal,T,max(sapply(xland,max))/4)
    mland <- meanLandscape(B,maxdim,maxscale,ret = T)
    m4 <- sapply(1:maxdim,function(d)count.local.maximal(mland[[d]],thresh = max(sapply(xland,max))/4))
    speak <- bootstrap.homology.mk2(B,maxdim,maxscale)
    m5 <- sapply(1:maxdim,function(d)speak[[paste0("dim",d,"dhole")]])
    
    aggr1[t,1] <- m1[1]
    aggr1[t,2] <- m2[1]
    aggr1[t,3] <- m3[1]
    aggr1[t,4] <- m4[1]
    aggr1[t,5] <- m5[1]
    
    aggr2[t,1] <- m1[2]
    aggr2[t,2] <- m2[2]
    aggr2[t,3] <- m3[2]
    aggr2[t,4] <- m4[2]
    aggr2[t,5] <- m5[2]
  }
  Xdiag<-lapply(1:length(X), function(k){return(X[[k]]$diag)})
  aggrs <- list(aggr1,aggr2)
  aggrs <- append(aggrs,list(Xsize=nrow(X[[1]]),Xsamples=length(X),
                             Bsize=size,Bsamples=samples,
                             maxdim=maxdim,maxscale=maxscale, Xdiag=Xdiag))
  class(aggrs) <- "bettiComp"
  
  return(aggrs)
}

#提案手法でのみサイクル数を推定する関数
proposedMethodOnly <- function(X,maxdim,maxscale,samples, const.size=0){
  aggr1 <- matrix(0,length(X),1)
  aggr2 <- matrix(0,length(X),1)
  dimnames(aggr1) <- list(paste0("data-set", 1:length(X)), "proposed")
  dimnames(aggr2) <- dimnames(aggr1)
  
  for(t in 1:length(X)){
    
    cat("data set", t, "calculating\n")
    if(const.size==0){size<-X[[t]]$nsample*(4/5)}
    else{size<-const.size}
    
    B <- bootstrapper(X[[t]]$noizyX,size,samples)
    # X[[t]]$diag<-calcPhom(X[[t]]$noizyX,maxdim,maxscale, plot = T,ret = T)
    # m1 <- calcDiagHomology(attr(B,"hausd"))[2:3]
    # m2 <- calcDiagHomology(hausdInterval(X[[t]]$noizyX,nrow(X[[t]]$noizyX)-1))[2:3]
    # xland <- lapply(1:maxdim,function(d)landscape(diagram[[1]],d))
    # m3 <- sapply(xland,count.local.maximal,T,max(sapply(xland,max))/4)
    # mland <- meanLandscape(B,maxdim,maxscale,ret = T)
    # m4 <- sapply(1:maxdim,function(d)count.local.maximal(mland[[d]],thresh = max(sapply(xland,max))/4))
    speak <- bootstrap.homology.mk2(B,maxdim,maxscale)
    m5 <- sapply(1:maxdim,function(d)speak[[paste0("dim",d,"dhole")]])
    
    # aggr1[t,1] <- m1[1]
    # aggr1[t,2] <- m2[1]
    # aggr1[t,3] <- m3[1]
    # aggr1[t,4] <- m4[1]
    # aggr1[t,5] <- m5[1]
    aggr1[t,1] <- m5[1]
    
    # aggr2[t,1] <- m1[2]
    # aggr2[t,2] <- m2[2]
    # aggr2[t,3] <- m3[2]
    # aggr2[t,4] <- m4[2]
    # aggr2[t,5] <- m5[2]
    aggr2[t,1] <- m5[2]
  }
  # Xdiag<-lapply(1:length(X), function(k){return(X[[k]]$diag)})
  aggrs <- list(aggr1,aggr2)
  # aggrs <- append(aggrs,list(Xsize=nrow(X[[1]]),Xsamples=length(X),
  #                            Bsize=size,Bsamples=samples,
  #                            maxdim=maxdim,maxscale=maxscale, Xdiag=Xdiag))
  aggrs <- append(aggrs,list(Xsize=sapply(1:length(X), function(l)nrow(X[[l]][["noizyX"]])),Xsamples=length(X),
                             Bsize=size,Bsamples=samples,
                             maxdim=maxdim,maxscale=maxscale))
  class(aggrs) <- "bettiComp"
  
  return(aggrs)
}

#提案手法でのみサイクル数を推定する関数。並列処理で
proposedMethodOnly.parallel <- function(X,maxdim,maxscale,samples){
  aggr1 <- matrix(0,length(X),1)
  aggr2 <- matrix(0,length(X),1)
  dimnames(aggr1) <- list(paste0("data-set", 1:length(X)), "proposed")
  dimnames(aggr2) <- dimnames(aggr1)
  
  for(t in 1:length(X)){
    
    cat("data set", t, "calculating\n")
    
    size<-X[[t]]$nsample*(3/5)
    B <- bootstrapper(X[[t]]$noizyX,size,samples)
    
    speak <- bootstrap.homology.parallel(B,maxdim,maxscale)
    m5 <- sapply(1:maxdim,function(d)speak[[paste0("dim",d,"dhole")]])
    
    aggr1[t,1] <- m5[1]
    
    aggr2[t,1] <- m5[2]
  }
  
  aggrs <- list(aggr1,aggr2)
  
  aggrs <- append(aggrs,list(Xsize=sapply(1:length(X), function(l)nrow(X[[l]][["noizyX"]])),Xsamples=length(X),
                             Bsize=size,Bsamples=samples,
                             maxdim=maxdim,maxscale=maxscale))
  class(aggrs) <- "bettiComp"
  
  
  
  return(aggrs)
}

#KDEをつかったPH図を計算
#KDE bootstrapを使った95%信頼区間によってサイクル判定
homMethodsKDE<- function(X, maxdim){
  require(TDA)
  aggr1 <- matrix(0,length(X),1)
  aggr2 <- matrix(0,length(X),1)
  dimnames(aggr1) <- list(paste0("data-set", 1:length(X)), "proposed")
  dimnames(aggr2) <- dimnames(aggr1)
  
  for(t in 1:length(X)){
    cat("data set", t, "calculating\n")
    
    Xlim <- range(X[[1]][["noizyX"]][,1])
    Ylim <- range(X[[1]][["noizyX"]][,2])
    Zlim <- range(X[[1]][["noizyX"]][,3])
    by <- (Xlim[2]-Xlim[1])/100
    Xseq <- seq(Xlim[1], Xlim[2], by = by)
    Yseq <- seq(Ylim[1], Ylim[2], by = by)
    Zseq <- seq(Zlim[1], Zlim[2], by = by)
    Grid <- expand.grid(Xseq, Yseq, Zseq)
    
    kde.diag<-gridDiag(X = X[[t]][["noizyX"]], FUN = kde, h = 0.3, lim = cbind(Xlim, Ylim, Zlim), by = by,
      sublevel = FALSE, library = "Dionysus", location = TRUE, printProgress = T)
    attr(kde.diag, "scale")<-c(0, max(data.collect1.DiagGrid[[1]][,3]))
    attr(kde.diag,"maxdimension")<-maxdim
    
    kde.band<-bootstrapBand(X = X[[t]][["noizyX"]], FUN = kde, h=0.3, Grid = Grid, B = 10, parallel = F, alpha = 0.05)
    m5 <- calcDiagHomology(kde.band[["width"]], kde.diag)[2:3]
    
    aggr1[t,1] <- m5[1]
    aggr2[t,1] <- m5[2]
    
    if(t==1){Diag<-kde.diag}
    else{Diag<-append(Diag, kde.diag)}
    
  }
  
  aggrs <- list(aggr1,aggr2)
  
  aggrs <- append(aggrs,list(Xsize=sapply(1:length(X), function(l)nrow(X[[l]][["noizyX"]])), Xsamples=length(X), maxdim=maxdim))
  class(aggrs) <- "bettiComp"

  return(aggrs)
}
