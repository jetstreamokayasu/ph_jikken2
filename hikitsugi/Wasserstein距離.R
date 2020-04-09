
calcWasserstein <- function(X,Y,maxdim,maxscale,Ymaxdim=maxdim,Ymaxscale=maxscale,wdim=1){
  Xdiag <- calcPhom(X,maxdim,maxscale,ret = T)$diagram
  Ydiag <- calcPhom(Y,Ymaxdim,Ymaxscale,ret = T)$diagram
  
  return(wasserstein(Xdiag,Ydiag,dimension = wdim))
}