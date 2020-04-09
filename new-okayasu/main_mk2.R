#二神先輩の関数であそぶためのスクリプト

#パーシステンスの平均を求める関数を変更
proposedMethodOnly.mk1 <- function(X,maxdim,maxscale,samples, const.size=0){
  aggr1 <- matrix(0,length(X),1)
  aggr2 <- matrix(0,length(X),1)
  dimnames(aggr1) <- list(paste0("data-set", 1:length(X)),"proposed")
  dimnames(aggr2) <- dimnames(aggr1)
  
  for(t in 1:length(X)){
    
    cat("data set", t, "calculating\n")
    if(const.size==0){size<-X[[t]]$nsample*(4/5)}
    B <- bootstrapper(X[[t]]$noizyX,size,samples)
    speak <- bootstrap.homology.mk3(B,maxdim,maxscale)
    m5 <- sapply(1:maxdim,function(d)speak[[paste0("dim",d,"dhole")]])
    
    #m6<-sapply(1:maxdim,function(d)mean(speak[["peak_am"]][d,]))
    #m7<-sapply(1:maxdim,function(d)mean(speak[["peak_zero"]][d,]))
    
    
    aggr1[t,1] <- m5[1]
    #aggr1[t,2] <- m6[1]
    #aggr1[t,3] <- m7[1]
    
    aggr2[t,1] <- m5[2]
    #aggr2[t,2] <- m6[2]
    #aggr2[t,3] <- m7[2]
  }
  # Xdiag<-lapply(1:length(X), function(k){return(X[[k]]$diag)})
  aggrs <- list(aggr1,aggr2)
  # aggrs <- append(aggrs,list(Xsize=nrow(X[[1]]),Xsamples=length(X),
  #                            Bsize=size,Bsamples=samples,
  #                            maxdim=maxdim,maxscale=maxscale, Xdiag=Xdiag))
  Xsize<-sapply(1:length(X), function(l){return(nrow(X[[l]][["noizyX"]]))})
  if(const.size==0){Bsize<-sapply(1:length(X), function(l){return(nrow(X[[l]][["noizyX"]])*(4/5))})}
  else{Bsize<-const.size}
  
  aggrs <- append(aggrs,list(Xsize=Xsize,Xsamples=length(X),
                             Bsize=Bsize,Bsamples=samples,
                             maxdim=maxdim,maxscale=maxscale))
  class(aggrs) <- "bettiComp"
  
  return(aggrs)
}

#パーシステンスの平均を求める関数を変更
bootstrap.homology.mk3 <- function(X,maxdim,maxscale,const.band=0,maximum.thresh = F){
  require(TDA)
  # require(pracma)
  if(!("bootsSamples" %in% class(X))) stop("input must be bootsSamples")
  peak <- matrix(0,maxdim,length(X))
  peak.am <- matrix(0,maxdim,length(X))
  #peak.zero <- matrix(0,maxdim,length(X))
  # band <- ifelse(const.band > 0,const.band,hausdInterval(X, m=sample.size, B=times, alpha = (1-confidence)))
  tseq <- seq(0,maxscale,length.out = 1000)
  diags <- lapply(X,function(x)calcPhom(x,maxdim,maxscale,ret = T,plot = F))
  print(sapply(diags,function(diag)calcPerMeanNoDouble(diag)))
  band <- ifelse(const.band==0,max(sapply(diags,function(diag)calcPerMeanNoDouble(diag))),const.band)
  print(band)
  
  #per.set<-lapply(diags, function(diag){return(calcPerdim1dim2(diag))})
  #band.am<-mean(unlist(per.set))
  #per0.set<-lapply(diags, function(diag){return(calcper(diag, 0))})
  #band.zero<-mean(unlist(per0.set))
    
    #debugText(band.am)
  
  for (t in 1:length(X)) {
    land <- lapply(1:maxdim,function(d)landscape(diags[[t]][[1]],dimension = d,KK = 1,tseq = tseq))
    if(maximum.thresh) band <- max(sapply(land,max))/4
    for(d in 1:maxdim){
      peak[d,t] <- calc.landscape.peak(X=land[[d]], thresh = (band/d), tseq=tseq)
      #peak.am[d,t] <- calc.landscape.peak(X=land[[d]], thresh = (band.am/d), tseq=tseq)
      #peak.zero[d,t] <- calc.landscape.peak(X=land[[d]], thresh = (band.zero/d), tseq=tseq)
    }
  }
  
  dimnames(peak) <- list(paste0("dim",1:maxdim),paste0("sample",1:length(X)))
  #dimnames(peak.am) <- list(paste0("dim",1:maxdim),paste0("sample",1:length(X)))
  #dimnames(peak.zero) <- list(paste0("dim",1:maxdim),paste0("sample",1:length(X)))
  bootstrap.summary <- list(peak=peak)
  bootstrap.summary <- append(bootstrap.summary,c(band=band,show.hole.density(peak)))
  class(bootstrap.summary) <- "smoothPhom"
  return(bootstrap.summary)
}