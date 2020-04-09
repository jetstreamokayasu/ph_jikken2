# mean landscape を描画します。
meanLandscape <- function(X,maxdim,maxscale,alpha = 0.95,col=c("red","pink","green","lightgreen"),ret=F){
  if(!have(class(X),"bootsSamples")) stop("argument X must be bootsSamples")
  diags <- lapply(X,function(x)CalcFastPhom(x,maxdim,maxscale,ret = T))
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
  class(X) <- c("matrix","bootsSamples")
  return(X)
}

# 4種の手法を比較
homologyMethodsComp <- function(X,maxdim,maxscale,size,samples){
  aggr1 <- matrix(0,5,20)
  aggr2 <- matrix(0,5,20)
  dimnames(aggr1) <- list(c("confidence","bestconf","landscape","meanland","proposed"),paste0(0:19,"holes"))
  dimnames(aggr2) <- dimnames(aggr1)
  
  for(t in 1:length(X)){
    B <- bootstrapper(X[[t]],size,samples)
    CalcFastPhom(X[[t]],maxdim,maxscale)
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
  aggrs <- append(aggrs,list(Xsize=nrow(X[[1]]),Xsamples=length(X),
                             Bsize=size,Bsamples=samples,
                             maxdim=maxdim,maxscale=maxscale))
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
    CalcFastPhom(X[[t]],maxdim,maxscale)
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