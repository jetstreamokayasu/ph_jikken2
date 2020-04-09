require(TDA)

rawdiag <- function(diag){
  if(class(diag)=="list") return(diag[["diagram"]])
  else return(diag)
}

calcPhom <- function(X, maxdimension, maxscale, plot = T, ret = F){
  diagram <<- ripsDiag(X,maxdimension = maxdimension,maxscale = maxscale,printProgress = T)
  if(plot) showPersistentDiagram(diagram)
  if(ret) return(diagram)
}

calcSubsamplePhom <- function(X, maxdimension, maxscale, rate = 0.5, plot = T, ret = F){
  assertthat::assert_that(assertthat::is.number(rate) && rate <= 1 && rate >= 0)
  subX <- X[sample(nrow(X),nrow(X)*rate),]
  calcPhom(subX, maxdimension, maxscale, plot, ret)
}

showPersistentDiagram <- function(diag = diagram,maxdimention = NA,scale = NA,point = 1.5){
  diag <- rawdiag(diag)
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

showPhomBand <- function(diag = diagram,band){
  diag <- rawdiag(diag)
  max.scale <- attr(diag,"scale")[2]
  
  plot(1,type = "n",ylim = c(0,max.scale),xlim = c(0,max.scale),ann = F,axes = F)
  polygon(c(max.scale+1,max.scale-band+1,-1,-1),c(max.scale+1,max.scale+1,band-1,0-1),col = "pink",border = 0)
  par(new=T)
  showPersistentDiagram(diag)
}

# パーシステンスが閾値以上のサイクルを表示
upperCycles <- function(band,diag = diagram,ignore.0 = F){
  diag <- rawdiag(diag)
  pdiag <- calcPersistence(diag)
  if(ignore.0) pdiag <- pdiag[pdiag[,"dimension"]!=0,]
  upperdiag <- pdiag[pdiag[,"Persistence"]>=band,]
  uppercycle <- upperdiag[order(upperdiag[,"Persistence"],decreasing = T),]
  if(nrow(uppercycle)==0) warning("there is no upper cycle")
  else return(uppercycle)
}

calcPersistence <- function(diag = diagram){
  diag <- rawdiag(diag)
  diag. <- cbind(diag,diag[,"Death"] - diag[,"Birth"])
  colnames(diag.) <- c(colnames(diag),"Persistence")
  return(diag.)
}
  
