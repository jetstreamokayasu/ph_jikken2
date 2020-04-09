library(pracma)

PhomFeature <- function(diag=diagmra,spar=0.7,weakcut=T) {
  diag <- diag[["diagram"]]
  tseq <- seq(0, attr(diag,"scale")[2], length.out = 1000)
  land <- landscape(diag,dimension = attr(diag,"maxdimension"),KK = 1,tseq = tseq)
  # if (missing(spar)){
  #   
  # }
  plot(tseq,land,type = "l",ylim = c(0,max(land)),col=2)
}

DiagramProcess <- function(X){
  if (class(X)!="list") {
    stop("diagram please")
  }
  d <- X[["diagram"]]
  dim0len <- length(which(d[,1]==0))
  out <- matrix(0,nrow(d),10)
  out[,1:3] <- d
  out[,4] <- Persistence(d)
  out[,5] <- c(rep(0,dim0len),Ranking(out[-(1:dim0len),4]))
  out[,6] <- c(rep(0,dim0len),Deviation(out[-(1:dim0len),4]))
  out[,7] <- sapply(X[["cycleLocation"]], length)
  out[,8] <- c(rep(0,dim0len),Deviation(out[-(1:dim0len),7]))
  out[,9] <- c(rep(0,dim0len),Deviation(sigmoid(out[-(1:dim0len),6]-max(out[,6]),a = 0.2)^2*out[-(1:dim0len),8]))
  out[,10]<- c(rep(0,dim0len),Ranking(out[-(1:dim0len),9]))
  colnames(out) <- c("dim","birth","death","persistence","prank","persisdev","cycles","cycledev","totaldev","totalrank")
  return(out)
}

Deviation <- function(X){
  z <- (X - mean(X)) / sd(X)
  s <- 50 + z * 10
  return(s)
}

Persistence <- function(X){
  if (class(X)=="list") {
    X <- X[["diagram"]]
  }
  return(X[,3]-X[,2])
}

Ranking <- function(X,col=NA){
  if (class(X)=="list") {
    X <- X[["diagram"]]
  }
  if (is.na(col)) {
    return(rank(-rank(X)))
  }
  else {
    return(rank(-rank(X[,col])))  
  }
}

NoisyExam <- function(X,smoothing = T){
  if (class(X)=="list") {
    X <- X[["diagram"]]
  }
  
  if (smoothing) {
    land <- smooth.spline(tseq,landscape(X,1,1,tseq),spar = 0.7)
    sil <- smooth.spline(tseq,silhouette(X,1,1,tseq),spar = 0.7)
  }
  else {
    land <- smooth.spline(tseq,landscape(X,1,1,tseq),spar = 0.01)
    sil <- smooth.spline(tseq,silhouette(X,1,1,tseq),spar = 0.01)
  }
  # print(paste("landscape spar =",as.character(land$spar)))
  # print(paste("silhouette spar =",as.character(sil$spar)))
  
  plot(land,col = 2,type = "l",xlab = "(birth+death)/2",ylab = "(death-birth)/2")
  par(new=T)
  plot(sil,col = 4,type = "l",ann = F,axes = F)
  abline(max(land$y)/4,0)
  
  print(DiffPeak(land$y))
  print(DiffPeak(sil$y))
  
}

MassTest <- function(n=100,spar = 0.7,savemode = F){
  landdata <- list()
  matdata <- matrix(0,n,2)
  for (i in 1:n) {
    pdiag <- ripsDiag(torusUnif(100,1,1),maxdimension = 1,maxscale = 2,library = 'Dionysus',printProgress = T)
    land <- landscape(pdiag[["diagram"]],dimension = 1,KK = 1,tseq = tseq)
    smooth <- smooth.spline(tseq,land,spar = spar)
    matdata[i,] <- c(DiffPeak(land),DiffPeak(smooth$y))
    
    plot(land,col = 2,type = "l",xlab = "(birth+death)/2",ylab = "(death-birth)/2")
    par(new=T)
    plot(smooth,col = 3,type = "l",ann = F,axes = F)
    abline(max(land)/4,0)
    if (savemode) {
      SavePNG(name = paste("land",as.character(i)," ",as.character(matdata[i,2]),"peaks",sep = ""))
    }
    
    print(paste(as.character(i),", smspeak = ",as.character(matdata[i,2]),sep = ""))
    
    landdata <- append(landdata,list(land))
  }
  return(matdata)
}

MassSumTest <- function(n=100,nthland = 1,spar = 0.7,savemode = F){
  sumlandscape <- rep(0,1000)
  matdata <- matrix(0,n,4)
  for (i in 1:n) {
    pdiag <- ripsDiag(torusUnif(100,1,1),maxdimension = 1,maxscale = 2,library = 'Dionysus',printProgress = T)
    land <- landscape(pdiag[["diagram"]],dimension = 1,KK = nthland,tseq = tseq)
    smooth <- smooth.spline(tseq,land,spar = spar)
    sumlandscape <- sumlandscape + land
    sumsmooth <- smooth.spline(tseq,sumlandscape,spar = spar)
    matdata[i,] <- c(DiffPeak(land),DiffPeak(smooth$y),DiffPeak(sumlandscape),DiffPeak(sumsmooth$y))
    
    plot(sumlandscape,col = 2,type = "l",xlab = "(birth+death)/2",ylab = "(death-birth)/2")
    par(new=T)
    plot(sumsmooth,col = 3,type = "l",ann = F,axes = F)
    abline(max(land)/4,0)
    if (savemode) {
      SavePNG(name = paste("land",as.character(i)," ",as.character(matdata[i,2]),"peaks",sep = ""))
    }
    
    print(paste(as.character(i),", smspeak = ",as.character(matdata[i,2]),sep = ""))
  }
  return(matdata)
}

# newdiag <- CalcPhom(torusUnif(100,1,1),maxdimention = 1,maxscale = 2)
# newland <- landscape(newdiag,dimension = 1,KK = 1,tseq = tseq)
# newsms <- 