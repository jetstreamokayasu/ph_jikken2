holeexistdiscri<-function(X){
  
  holeexist<-matrix(0, length(X), 2)
  
  for(i in 1:length(X)){
    
    per.dim0<-X[[i]]$diag[X[[i]]$diag[,1]==0, 3]
    per.dim0<-per.dim0[-1]
    mean.dim0<-mean(per.dim0)
    
    per.dim1<-X[[i]]$diag[X[[i]]$diag[,1]==1, 3]-X[[i]]$diag[X[[i]]$diag[,1]==1, 2]
    
    per.dim2<-X[[i]]$diag[X[[i]]$diag[,1]==2, 3]-X[[i]]$diag[X[[i]]$diag[,1]==2, 2]
    per.dim2<-per.dim2*2
    
    cat("data set", i, "\n")
    cat("mean.dim0=", mean.dim0, "\n")
    cat("max.per.dim1=", per.dim1[which.max(per.dim1)], "\n")
    cat("max.per.dim2=", per.dim2[which.max(per.dim2)], "\n")
    
    
    
    exist.dim1<-per.dim1[per.dim1>=mean.dim0]
    exist.dim2<-per.dim1[per.dim1>=mean.dim0]
    
    if(length(exist.dim1)>0){
      
      holeexist[i, 1]<-1
      
      cat("data", i, "dim1 hole exist \n")
      
    }
    
    else{cat("data", i, "no dim1 hole exist \n")}
    
    
    if(length(exist.dim2)>0){
      
      holeexist[i, 2]<-1
      
      cat("data", i, "dim2 hole exist \n")
      
    }
    else{cat("data", i, "no dim2 hole exist \n")}
    
    
    
  }
  
  return(holeexist)
  
  
}


holexistdetermine<-function(X){
  #サイクル数推定結果のPH図からサイクルがあるか判定
  #閾値は0次パーシステンスの平均
  
  holeexist<-matrix(0, length(X), 2)
  
  for(i in 1:length(X)){
    
    per.dim0<-X[[i]]$diagram[X[[i]]$diagram[,1]==0, 3]
    per.dim0<-per.dim0[-1]
    mean.dim0<-mean(per.dim0)
    
    #debugText(mean.dim0)
    
    per.dim1<-X[[i]]$diagram[X[[i]]$diagram[,1]==1, 3]-X[[i]]$diagram[X[[i]]$diagram[,1]==1, 2]
    
    per.dim2<-X[[i]]$diagram[X[[i]]$diagram[,1]==2, 3]-X[[i]]$diagram[X[[i]]$diagram[,1]==2, 2]
    per.dim2<-per.dim2*2
    
    cat("data set", i, "\n")
    cat("mean.dim0=", mean.dim0, "\n")
    cat("max.per.dim1=", per.dim1[which.max(per.dim1)], "\n")
    cat("max.per.dim2=", per.dim2[which.max(per.dim2)], "\n")
    cat("mean.per.dim1=", mean(per.dim1), "\n")
    cat("mean.per.dim2=", mean(per.dim2), "\n")
    cat("mean.per=", mean(c(per.dim1,per.dim2)), "\n")
    cat("mean.per*2=", mean(c(per.dim1,per.dim2))*2, "\n")
    
    
    exist.dim1<-per.dim1[per.dim1>=mean.dim0]
    exist.dim2<-per.dim1[per.dim2>=mean.dim0]
    
    if(length(exist.dim1)>0){
      
      holeexist[i, 1]<-1
      
      cat("data", i, "dim1 hole exist \n")
      
    }
    
    else{cat("data", i, "no dim1 hole exist \n")}
    
    
    if(length(exist.dim2)>0){
      
      holeexist[i, 2]<-1
      
      cat("data", i, "dim2 hole exist \n")
      
    }
    else{cat("data", i, "no dim2 hole exist \n")}
    
    
    
  }
  
  cat(length(which(holeexist[,1]==1)), "sets of data have a 1-dim hole\n")
  cat(length(which(holeexist[,2]==1)), "sets of data have a 2-dim hole\n")
  
  cycle.set<-c(length(which(holeexist[,1]==1)), length(which(holeexist[,2]==1)))
  names(cycle.set)<-c("1-dim", "2-dim")
  
  return(list(holeexist=holeexist, cycle.set=cycle.set))
  
  
}

#サイクル数がn個あると判されたデータセットがそれぞれいくつあるか調べる
#信頼区間、PL、MLと比較する場合とそれ以外で場合分け
cyclenumber<-function(holes, compare=T, method=0){
  
  cyc.num<-max(holes)
  #debugText(compare)
  
  if(cyc.num > 4){col<-cyc.num}
  else{col<-4}
  
  if(compare==T){
  estimate<-matrix(0, 5, col+1)
  rownames(estimate)<-c("confidence", "bestconf", "landscape", "meanland", "proposed")
  colnames(estimate)<-c(paste0(0:col, "hole"))
  
  #debugText(estimate)
  
  for (i in 1:(col+1)) {
    #debugText(i)
    estimate["confidence", i]<-meetNumber(holes[, "confidence"], i-1)
    estimate["bestconf", i]<-meetNumber(holes[, "bestconf"], i-1)
    estimate["landscape", i]<-meetNumber(holes[, "landscape"], i-1)
    estimate["meanland", i]<-meetNumber(holes[, "meanland"], i-1)
    estimate["proposed", i]<-meetNumber(holes[, "proposed"], i-1)
    
  }
  
  }
  else{
      estimate<-matrix(0, ncol(holes), col+1)
      if(method[1]==0){rownames(estimate)<-c(paste0("method", 1:ncol(holes)))}
      else{rownames(estimate)<-method}
      colnames(estimate)<-c(paste0(0:col, "hole"))
      
      #debugText(estimate)
      
      for (i in 1:(col+1)) {
        for (j in 1:ncol(holes)) {
          estimate[j, i]<-meetNumber(holes[, j], i-1)
        }
      }
    
  }
  
  return(estimate)
  
}

meetNumber<-function(result, correct){
  
  return(length(result[result >= (correct-0.5) & result < (correct+0.5)]))
  
}

cyclenumber2<-function(holes, exist, dim){
  #サイクルが存在すると判定されたデータセットのサイクル推定数のみ計算
  
  estimate<-matrix(0, 5, 5)
  rownames(estimate)<-c("confidence", "bestconf", "landscape", "meanland", "proposed")
  colnames(estimate)<-c("0hole", "1hole", "2hole", "3hole", "4hole")
  
  for (i in 1:5) {
    
    estimate["confidence", i]<-meetNumber(holes[which(exist[, dim]==1), "confidence"], i-1)
    estimate["bestconf", i]<-meetNumber(holes[which(exist[, dim]==1), "bestconf"], i-1)
    estimate["landscape", i]<-meetNumber(holes[which(exist[, dim]==1), "landscape"], i-1)
    estimate["meanland", i]<-meetNumber(holes[which(exist[, dim]==1), "meanland"], i-1)
    estimate["proposed", i]<-meetNumber(holes[which(exist[, dim]==1), "proposed"], i-1)
    
  }
  
  return(estimate)
  
}

holexistdetermine2<-function(X){
  #サイクル数推定結果のPH図からサイクルがあるか判定
  #閾値は指数分布に基づいて設定
  
  if(class(X)=="diagram"){X<-list(X)}
  
  holeexist<-matrix(0, length(X), 2)
  
  for(i in 1:length(X)){
    
    per.dim1<-calcper(X[[i]],1)
    per.dim2<-calcper(X[[i]],2)
    per.dim2<-per.dim2*2
    
    thresh<-threshDetermine(X[[i]], 3)
    
    cat("data set", i, "\n")
    #cat("mean.dim0=", mean.dim0, "\n")
    cat("max.per.dim1=", per.dim1[which.max(per.dim1)], "\n")
    cat("max.per.dim2=", per.dim2[which.max(per.dim2)], "\n")
    #cat("mean.per.dim1=", mean(per.dim1), "\n")
    #cat("mean.per.dim2=", mean(per.dim2), "\n")
    #cat("mean.per=", mean(c(per.dim1,per.dim2)), "\n")
    #cat("mean.per*2=", mean(c(per.dim1,per.dim2))*2, "\n")
    cat("thresh=", thresh, "\n")
    
    exist.dim1<-per.dim1[per.dim1>=thresh]
    exist.dim2<-per.dim1[per.dim2>=thresh]
    
    if(length(exist.dim1)>0){
      
      holeexist[i, 1]<-1
      
      cat("data", i, "dim1 hole exist \n")
      
    }
    
    else{cat("data", i, "no dim1 hole exist \n")}
    
    
    if(length(exist.dim2)>0){
      
      holeexist[i, 2]<-1
      
      cat("data", i, "dim2 hole exist \n")
      
    }
    else{cat("data", i, "no dim2 hole exist \n")}
    
    
    
  }
  
  cat(length(which(holeexist[,1]==1)), "sets of data have a 1-dim hole\n")
  cat(length(which(holeexist[,2]==1)), "sets of data have a 2-dim hole\n")
  
  cycle.set<-c(length(which(holeexist[,1]==1)), length(which(holeexist[,2]==1)))
  names(cycle.set)<-c("1-dim", "2-dim")
  
  return(list(holeexist=holeexist, cycle.set=cycle.set))
  
  
}

calcper<-function(diag, dim){
  #dimに指定した次元のそれぞれの穴のパーシステンスを求める関数
  if(class(diag)=="list"){diag<-diag[[1]]}
  per.dim<-diag[diag[,1]==dim, 3]-diag[diag[,1]==dim, 2]
  #rownames(per.dim)<-paste0(dim, "per")
  if(dim==0){per.dim<-per.dim[-1]}
  
  return(per.dim)
}


threshDetermine<-function(diag, maxscale){
  #サイクルの有る無しを判定する閾値を決める関数
  #指数分布に基づいて設定
  diag<-rawdiag(diag)
  #maxdim<-max(diag[,1])
  
  persis<-lapply(1:2, function(l)calcper(diag, l))
  
  allpersis<-c(persis[[1]], persis[[2]]*2)
  
  #debugText(length(allpersis))
  
  thresh<-qexp(0.99, (length(allpersis)/maxscale))*maxscale
  
  return(thresh)
  
}

#1次元と2次元のパーシステンスを計算、1つのベクトルにする
calcPerdim1dim2<-function(diag){
  
  if(class(diag)=="list"){diag<-diag[[1]]}
  per.dim1<-calcper(diag, 1)
  per.dim2<-calcper(diag, 2)
  
  per.dim1dim2<-c(per.dim1, per.dim2*2)
  
  return(per.dim1dim2)
  
}
