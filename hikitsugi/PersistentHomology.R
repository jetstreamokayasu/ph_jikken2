######################################################
##TDA Support
######################################################
library(TDA)

showTorus <- function(X){
  require(rgl)
  plot3d(X, aspect=F,xlab = "",ylab = "",zlab = "",cex = 1.2)
}

calcLand <- function(X, maxdimention = 1, KK = 1){
  if (class(X)=="list") X <- X[["diagram"]]
  if (exists("tseq")==F) tseq <- seq(0,attr(X,"scale")[2],length.out = 1000)
  #if (missing(maxdimention)) maxdimention <- attr(diag,"maxdimension")
  
  land <- landscape(X,dimension = maxdimention,KK = KK,tseq = tseq)
  return(land)
}

calcDiagCentroid <- function(diag = diagram,target.dim = 1){
  if(class(diag)=="list") diag <- diag[[1]]
  maxdim <- attr(diag,"maxdimension")
  ind <- setdiff(1:nrow(diag),which(diag[,1]==setdiff(0:maxdim,target.dim)))
  diag <- diag[ind,]
  
  centroid <- apply(diag[,2:3],2,mean)
  cpersistence <- (centroid[2]-centroid[1])
  ret <- c(centroid,cpersistence,cpersistence*2)
  if(nrow(diag)==0) ret <- numeric(4)
  names(ret) <- c("birth","death","persistence","noizes.thresh")
  return(ret)
}

calcFastPhom <- function(X, maxdimension, maxscale,plot = T,ret = F){
  diagram <<- ripsDiag(X,maxdimension = maxdimension,maxscale = maxscale,printProgress = T)
  if(plot) showPHOM(diagram)
  if(ret) return(diagram)
}

showPHOM <- function(diag=diagram,maxdimention=NA,scale=NA,point=1.5){
  if(class(diag)=="list") diag <- diag[["diagram"]]
  if(missing(maxdimention)) maxdimention <- attr(diag,"maxdimension")
  if(missing(scale)) scale <- attr(diag,"scale")
  maxdimention <- maxdimention + 1
  diag[is.infinite(diag)] <- scale[2]
  
  plot(diag[,2], diag[,3], xlim=scale, ylim=scale,	
       cex=point, cex.axis=point, col=diag[,1] + 1,
       pch=diag[,1] + 1, xlab="Birth", ylab="Death",
       cex.lab=point, cex.main=2)
  abline(0,1)
  
  legends <- 0
  for (i in (0:(maxdimention-1))){
    # legends <- c(legends,bquote(H[.(i)])) 
    legends <- c(legends,paste("dim",i))
  }
  legends <- legends[-1]
  
  legend(scale[2]/2*1.2,scale[2]/2,legend=sapply(legends,as.expression),
        col=1:maxdimention,pch=1:maxdimention,cex=point,pt.cex=point)
}

showPhomWithBand <- function(band,diag=diagram){
  diag.mat <- diag[[1]]
  max.scale <- attr(diag.mat,"scale")[2]
  plot(1,type = "n",ylim = c(0,max.scale),xlim = c(0,max.scale),ann = F,axes = F)
  polygon(c(max.scale+1,max.scale-band+1,-1,-1),c(max.scale+1,max.scale+1,band-1,0-1),col = "pink",border = 0)
  par(new=T)
  showPHOM(diag = diag)
}

showLandscape <- function(diag=diagram,dimension = 1, KK = 1,ylimit = -1){
  if(exists("tseq")==F) tseq <- seq(0,attr(diag[["diagram"]],"scale")[2],length.out = 1000)
  land <- landscape(diag[["diagram"]], dimension = dimension, KK = KK, tseq = tseq)
  if (ylimit == -1){
    ylim <- c(0,max(land))
  }
  else {
   ylim <- c(0,ylimit) 
  }
  plot(tseq,land,type = "l",ylim = ylim)
}

showSilhouette <- function(diag = diagram, dimension = 1, p = 1,ylimit = -1){
  sil <- silhouette(diag[["diagram"]], dimension = dimension, p = p, tseq = tseq)
  if (ylimit == -1){
    ylim <- c(0,max(sil))
  }
  else {
    ylim <- c(0,ylimit) 
  }
  plot(tseq,sil,type = "l",ylim = ylim)
}


showSumLandscape <- function(diag = diagram, maxdimention = 1, KK = 1,cut = F){
  if(exists("tseq")==F) tseq <- seq(0,attr(diag[["diagram"]],"scale")[2],length.out = 1000)
  if(missing(maxdimention)) maxdimention <- attr(diag[["diagram"]],"maxdimension")
  lands <- lapply(1:maxdimention,function(d)landscape(diag[["diagram"]], dimension = d, KK = KK, tseq = tseq))
  lim <- c(0,max(sapply(1:maxdimention,function(d)max(lands[[d]]))))
  for (dim in 1:maxdimention){
    if (dim==1) plot(tseq,lands[[dim]],type = "l",col=dim+1,ylim = lim,xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2")
    else plot(tseq,lands[[dim]],type = "l",col=dim+1,ylim = lim,ann = F,axes = F)
    par(new=T)
  }
  
  legend(0,lim[2],legend=sapply(1:maxdimention,function(d)paste("dim",d)),
         col=2:(maxdimention+1),pch=3)
  
  if(cut > 0) abline(lim[2]/4,0)
  par(new=F)
  return(lim)
}