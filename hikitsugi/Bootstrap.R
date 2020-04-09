GenerateSample <- function(x,samples=100){
  if (!is.matrix(x)) stop("input matrix data")
  return(x[order(runif(nrow(x)))[1:samples],])
}

BootstrapPhom <- function(x,maxdim,maxscale,nsample=100,ntime=10,spar=seq(0,1,0.1)){
  peakstats <- list()
  for (i in 1:ntime) {
    CalcFastPhom(GenerateSample(x,nsample),maxdim,maxscale)
    land <- CalcLand(diagram)
    tseq <- seq(0,maxscale,length.out = 1000)
    peaks <- list(DiffPeak(land,tseq = tseq))
    peaks <- c(peaks,lapply(spar, function(sp){SmoothDiffPeak(land,sp,tseq=tseq)}))
    names(peaks) <- c("normal",as.character(spar))
    peakstats <- append(peakstats,list(peaks))
  }
  peakstats <<- peakstats
}

showBootstrap <- function(stats = peakstats){
  x <- matrix(0,nrow = length(stats),ncol = length(stats[[1]]))
  dataset <- lapply(peakstats, function(peaks){return(unlist(peaks))})
  
  for (i in 1:length(stats)) {
    for (j in 1:length(stats[[1]])) {
      x[i,j] <- x[i,j] + stats[[i]][[j]]
    }
  }
  colnames(x) <- names(stats[[1]])
  barplot(x)
  return(x)
}

# ShowFeatureProbability <- function(X = showBootstrap(peakstats)) {
#   maxpeak <- max(X)
#   prob <- numeric(maxpeak+1)
#   for (i in 0:maxpeak) {
#     prob[i] <- prob[i] + length(which(X==i))
#   }
#   return(prob/lengtaah(X))
# }

#平滑化してピーク数を数える。ダイアグラム直入力可。
SmoothDiffPeak <- function(X = diagram,spar = 0.7,tseq = tseq){
  if (class(X)=="list") {
    X <- X[["diagram"]]
    tseq  <- seq(0,attr(X,"scale")[2],length.out = 1000)
    X <- landscape(X,1,1,tseq)
  }
  SMS <- smooth.spline(tseq,X,spar = spar)
  plot(SMS,type = "l")
  abline(max(SMS$y)/4,0)
  return(DiffPeak(SMS$y,tseq = tseq))
}

#landscapeのピーク数を数える
DiffPeak <- function(vec,weakcut = T,thresh = max(vec)/4,tseq = tseq){
  peak <- 0
  d <- diff(vec)
  thresh <- thresh * as.numeric(weakcut)
  for (i in 1:(length(vec)-2)) {
    if (d[i]*d[i+1] < 0 && vec[i] > thresh && d[i]>d[i+1]) {
      peak <- peak + 1
    }
  }
  return(peak)
}