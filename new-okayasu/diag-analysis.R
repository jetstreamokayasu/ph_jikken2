hutagamiHomologyMethodHutagami <- function(X,maxdim,maxscale,size.set,samples){
  aggr1 <- matrix(0,5,20)
  aggr2 <- matrix(0,5,20)
  dimnames(aggr1) <- list(c("confidence","bestconf","landscape","meanland","proposed"),paste0(0:19,"holes"))
  dimnames(aggr2) <- dimnames(aggr1)
  
  for(t in 1:length(X)){
    size<-size.set[t]*(3/5)
    B <- bootstrapper(X[[t]],size,samples)
    calcPhom(X[[t]],maxdim,maxscale)
    m1 <- calcDiagHomology(attr(B,"hausd"))[2:3]
    m2 <- calcDiagHomology(hausdInterval(X[[t]],nrow(X[[t]])-1))[2:3]
    xland <- lapply(1:maxdim,function(d)landscape(diagram[[1]],d))
    m3 <- sapply(xland,count.local.maximal,T,max(sapply(xland,max))/4)
    mland <- meanLandscape(B,maxdim,maxscale,ret = T)
    m4 <- sapply(1:maxdim,function(d)count.local.maximal(mland[[d]],thresh = max(sapply(xland,max))/4))
    speak <- bootstrap.homology(B,maxdim,maxscale)
    m5 <- sapply(1:maxdim,function(d)round(speak[[paste0("dim",d,"dhole")]]))
    cat("m5[1]=", m5[1], "\n")
    aggr1[1,m1[1]+1] <- aggr1[1,m1[1]+1] + 1
    aggr1[2,m2[1]+1] <- aggr1[2,m2[1]+1] + 1
    aggr1[3,m3[1]+1] <- aggr1[3,m3[1]+1] + 1
    aggr1[4,m4[1]+1] <- aggr1[4,m4[1]+1] + 1
    aggr1[5,m5[1]+1] <- aggr1[5,m5[1]+1] + 1
    hole.set[t, 3]<<-m5[1]
    
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