

# normalizing<-function(data, max, min){
#   
#   max.data<-max(data)
#   min.data<-min(data)
#   
#   nized.data<-sapply(data, function(x){return((((x-min.data)/(max.data-min.data))*(max-min))+min)})
#   
#   return(nized.data)
#   
# }

uniDisMake<-function(N, R){
  #半径Rの円内に一様分布する点群を作る関数
  theta <- runif(N, min=0, max=2*pi)
  r <- sqrt(2*runif(N, min=0, max=0.5*R^2))
  df <- data.frame(x=r*cos(theta), y=r*sin(theta))
  #ggplot(df, aes(x=x, y=y)) + geom_point(color="blue")
  return(df)
  
  
}

#ランドスケープを計算
calcLandscape<-function(diag, line=T){
  
  maxscale <- 3
  thresh<-calcDiagCentroid.mk2(diag)[1]
  thresh.kde<-threshPerKDE(diag, 0.9)
  thresh.zero<-mean(calcper(diag, 0))
  tseq <- seq(0, maxscale, length = 1000) #domain
  Land.dim1 <- landscape(diag[[1]], dimension = 1, KK = 1, tseq)
  plot(tseq, Land.dim1, type = "l", col=2, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, round(max(Land.dim1)+1)/2), main ="1-degree landscape")
  if(line){
  abline(h=thresh)
  #abline(h=thresh.kde/2, col=4)
  #abline(h=thresh.zero, col="orange")
  }
  
  if(length(diag[[1]][diag[[1]][,1]==2,])>0){
    Land.dim2 <- landscape(diag[[1]], dimension = 2, KK = 1, tseq)
    #if(new){par(new=T)}
    plot(tseq, Land.dim2, type = "l", col=3, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, round(max(Land.dim2)+1)/2), main ="2-degree landscape")
    if(line){
    abline(h=thresh/2)
    #abline(h=thresh.kde/4, col=4)
    #abline(h=thresh.zero/2, col="orange")
    }
    
    return(list(tseq=tseq, Land.dim1=Land.dim1, Land.dim2=Land.dim2, thresh=thresh))
    
  }else{return(list(tseq=tseq, Land.dim1=Land.dim1, thresh=thresh))}
}

#ランドスケープを描写
plotLandscape<-function(land){
  
  plotland<-lapply(2:(length(land)-1), function(k){
    
    if(names(land)[k]=="Land.dim1"){
    
    plot(land[[1]], land[[k]], type = "l", col=k, xlab = "(Birth + Death) / 2", ylab = "(Death - Birth) / 2", ylim=c(0, round(max(land[[k]])+1)/2), main =paste0(k-1, "-degree landscape"))
    abline(h=land[["thresh"]])
      
    }
    
    else if(names(land)[k]=="Land.dim2"){
      
      plot(land[[1]], land[[k]], type = "l", col=3, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, round(max(land[[k]])+1)/2), main =paste0(2, "-degree landscape"))
      abline(h=land[["thresh"]]/2)
      
    }else{
      
      plot(land[[1]], land[[k]], type = "l", col=k, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, round(max(land[[k]])+1)/2), main =paste0(k-1, "-degree landscape"))
      
    }
    
    })
  
}

#ワイブル分布による閾値設定
threshWeibull<-function(diag){
  
  diag<-rawdiag(diag)
  dim1.per<-calcper(diag, 1)
  dim2.per<-calcper(diag, 2)
  all.per<-c(dim1.per, dim2.per*2)
  debugText(length(dim1.per), length(dim2.per), length(all.per))
  wei<-weibull.par(all.per, epsilon=1e-7)
  thresh<-qweibull(0.95, shape = wei[1], scale = wei[2])
  
  return(thresh)
  
}

#3次元図形をプロット
figurePlot<-function(X){
  
  require(rgl)
  plot3d(X)
  aspect3d("iso")
  
}

#TDAパッケージの例から持ってきた関数
supnor.bottle <- function(Diag, data, indices){
  #debugText(indices)
  #indices<-sample(nrow(data), size=(nrow(data)*(3/5)))
  data.extract <- data[indices,]
  kde.Diag.boot <- gridDiag(X = data.extract, FUN = kde, h = 0.3,
                                   lim = cbind(range(data.extract[,1]), range(data.extract[,2]), range(data.extract[,3])),by = 0.065,
                                   sublevel = FALSE, library = "Dionysus", location = TRUE,
                                   printProgress = FALSE)
  
  bottle <- bottleneck(kde.Diag.boot[["diagram"]][kde.Diag.boot[["diagram"]][,1]==1,], Diag)
  return(bottle)
}

supnor <- function(kde, data,indices){
  data.extract <- data[indices,]
  by<-(range(data.extract[,1])[2]-range(data.extract[,1])[1])/100
  Xseq <- seq(range(data.extract[,1])[1], range(data.extract[,1])[2], by = by)
  Yseq <- seq(range(data.extract[,2])[1], range(data.extract[,2])[2], by = by)
  Zseq <- seq(range(data.extract[,3])[1], range(data.extract[,3])[2], by = by)
  Grid <- expand.grid(Xseq, Yseq, Zseq)
  hatf <- kde(X = data.extract, Grid = Grid, h = 0.3)
  supnorm <- max(abs(hatf - kde))
  return(supnorm)
}

#パーシステンスの分布をKDEによって推定して95%ラインを閾値に
threshPerKDE<-function(diag, line){
  
  require(ks)
  require(MaskJointDensity)
  
  diag<-rawdiag(diag)
  
  per.dim1dim2<-calcPerdim1dim2(diag)
  length(per.dim1dim2)

  per.kde<-ks::kde(sort(per.dim1dim2))
  thre<-qkdeSorted(line, per.kde)

  return(thre)
  
}

# threshPerKDE.norequire<-function(diag){
#   
#   diag<-rawdiag(diag)
#   
#   per.dim1dim2<-calcPerdim1dim2(diag)
#   
#   per.kde<-ks::kde(sort(per.dim1dim2))
#   thre<-qkdeSorted(p = 0.95, fhat = per.kde)
#   
#   return(thre)
#   
# }


#二神先輩の関数で遊ぶ
#提案手法でのみサイクル数を推定する関数 with KDE閾値
proposedMethodOnly.kde <- function(X,maxdim,maxscale,samples, const.size=0){
  #if("package:TDA" %in% search()){detach(package:TDA)}
  require(ks)
  require(MaskJointDensity)
  aggr1 <- matrix(0,length(X),1)
  aggr2 <- matrix(0,length(X),1)
  dimnames(aggr1) <- list(paste0("data-set", 1:length(X)), "proposed")
  dimnames(aggr2) <- dimnames(aggr1)
  
  for(t in 1:length(X)){
    
    cat("data set", t, "calculating\n")
    if(const.size==0){size<-X[[t]]$nsample*(4/5)}
    B <- bootstrapper(X[[t]]$noizyX,size,samples)
    X[[t]]$diag<-calcPhom(X[[t]]$noizyX,maxdim,maxscale, plot = T,ret = T)
    thresh<-threshPerKDE(X[[t]]$diag, 0.9)
    speak <- bootstrap.homology.kde(B,maxdim,maxscale, const.band = thresh)
    m5 <- sapply(1:maxdim,function(d)speak[[paste0("dim",d,"dhole")]])
    
    aggr1[t,1] <- m5[1]
   
    aggr2[t,1] <- m5[2]
  }
  Xdiag<-lapply(1:length(X), function(k){return(X[[k]]$diag)})
  aggrs <- list(aggr1,aggr2)
  aggrs <- append(aggrs,list(Xsize=sapply(1:length(X), function(l)nrow(X[[l]][["noizyX"]])),Xsamples=length(X),
                             Bsize=size,Bsamples=samples,
                             maxdim=maxdim,maxscale=maxscale, Xdiag=Xdiag))
  class(aggrs) <- "bettiComp"
  
  return(aggrs)
}

#KDEによる閾値設定を試す
bootstrap.homology.kde <- function(X,maxdim,maxscale, const.band=0, maximum.thresh = F){
  require(TDA)
  # require(pracma)
  if(!("bootsSamples" %in% class(X))) stop("input must be bootsSamples")
  peak <- matrix(0,maxdim,length(X))
  # band <- ifelse(const.band > 0,const.band,hausdInterval(X, m=sample.size, B=times, alpha = (1-confidence)))
  tseq <- seq(0,maxscale,length.out = 1000)
  diags <- lapply(X,function(x)calcPhom(x,maxdim,maxscale,ret = T,plot = F))
  #print(sapply(diags,function(diag)threshPerKDE.norequire(diag)))
  band <- ifelse(const.band==0,max(sapply(diags,function(diag)threshPerKDE.norequire(diag))),const.band)
  print(band)
  
  
  for (t in 1:length(X)) {
    land <- lapply(1:maxdim,function(d)landscape(diags[[t]][[1]],dimension = d,KK = 1,tseq = tseq))
    if(maximum.thresh) band <- max(sapply(land,max))/4
    for(d in 1:maxdim){
      peak[d,t] <- calc.landscape.peak(X=land[[d]], thresh = band/(2*d), tseq=tseq)
    }
  }
  
  dimnames(peak) <- list(paste0("dim",1:maxdim),paste0("sample",1:length(X)))
  bootstrap.summary <- list(peak=peak)
  bootstrap.summary <- append(bootstrap.summary,c(band=band,show.hole.density(peak)))
  class(bootstrap.summary) <- "smoothPhom"
  return(bootstrap.summary)
}

#データコレクトからサブサンプルセットを作る
subsampleExclude<-function(collect, nsub){
  
  new.collect<-lapply(collect, function(set){
    data<-set[["noizyX"]][sample(nrow(set[["noizyX"]]), nsub),]
    return(list(nsample=nsub, noizyX=data, diag=0))
  })
  
  return(new.collect)
  
}

#データセットcollectのそれぞれのデータセットに対し
#データ点nsub個のサブサンプルを作成
#それに対してベッチ数推定を行う
#これをtimes回繰り返す
successTimes<-function(collect, nsub, times, maxdim, maxscale, subsample=T){
  
  aggr.list<-lapply(1:times, function(t){
    
    cat(t, "times subcollect\n")
    
    if(subsample){
    subcollect<-subsampleExclude(collect, nsub)
    aggr<-proposedMethodOnly.mk1(subcollect, maxdim, maxscale, 10)
    }else{
      aggr<-proposedMethodOnly.mk1(collect, maxdim, maxscale, 10)
    }
    
    cycle<-lapply(1:maxdim, function(dim){
      
      result<-cyclenumber(aggr[[dim]], compare = F)
      cat(t, "times subcollect\n")
      cat(dim, "dimension\n")
      print(result)
      return(result)
      
      })
    
    names(cycle)<-paste0("dim", 1:maxdim, "cycles")
    #debugText(aggr)
    return(append(aggr, cycle))
    
    
  })
  
  return(aggr.list)
  
}

#パーシステンスの平均を求める
calcPerMean<-function(diag){
  
  if(class(diag)=="list"){diag<-diag[[1]]}
  
  if(2 %in% diag[,1]){
    per.mean<-mean(calcPerdim1dim2(diag))
  }else{
    per.mean<-mean(calcper(diag, 1))
  }
  
  names(per.mean)<-"mean_per"
  
  return(per.mean)
  
}

#パーシステンスの平均を求める
#2次パーシステンスを2倍にしない
calcPerMeanNoDouble<-function(diag){
  
  if(class(diag)=="list"){diag<-diag[[1]]}
  
  if(2 %in% diag[,1]){
    per.dim1<-calcper(diag, 1)
    per.dim2<-calcper(diag, 2)
    per.mean<-mean(c(per.dim1, per.dim2))
    
  }else{
    per.mean<-mean(calcper(diag, 1))
  }
  
  names(per.mean)<-"mean_per"
  
  return(per.mean)
  
}

#それぞれのサブサンプルセットごとの成功率を見る
successRates<-function(aggr.list, dim, correct, inter=F){
  
  if(inter==F){rates<-lapply(aggr.list, function(aggr){
    
    set.sum<-sum(aggr[[paste0("dim", dim, "cycles")]])
    
    rate<-aggr[[paste0("dim", dim, "cycles")]][correct+1]/set.sum
    
    return(rate)
    
  })}
  
  else{rates<-lapply(aggr.list, function(aggr){
    set.sum<-sum(aggr[[dim+18]])
    rate<-aggr[[dim+18]][correct+1]/set.sum
    return(rate)
  })}
  
  if(inter==F){names(rates)<-c(paste0("rate", 1:length(aggr.list)))}
  else{names(rates)<-c(paste0("inrate", 1:length(aggr.list)))}
  
  return(rates)
  
}

#aggrから直接成功率を求める
aggrSuccessRates<-function(aggr.list, correct){
  
  rates<-lapply(aggr.list, function(aggr){
    
    set.sum<-length(aggr[[1]])
    
    cycle1<-cyclenumber(aggr[[1]], compare=F)
    rate1<-cycle1[correct[1]+1]/set.sum
    
    cycle2<-cyclenumber(aggr[[2]], compare=F)
    rate2<-cycle2[correct[2]+1]/set.sum
    
    return(list(dim1rate=rate1, dim2rate=rate2))
    
  })

  
  return(rates)
  
}

#collectデータセットのについて段階的にサブサンプルを取りながら成功率を見る関数
#collece=オリジナルデータ、nsubset=サブサンプル数のセット、times=サブサンプルを何回行うか
succesCheck<-function(collect, nsubset, times, maxdim, mxscale, corrects, origin=F){
  
  rates.list<-lapply(nsubset, function(nsub){
    cat(nsub, "subsample calculating\n")
    aggr.list<-successTimes(collect, nsub, times, maxdim, mxscale)
    rates<-lapply(1:maxdim, function(dim){
      return(successRates(aggr.list, dim, corrects[dim]))
    })
    cat(nsub, "subsample\n")
    names(rates)<-c(paste0("dim", 1:maxdim))
    print(rates)
    return(rates)
  })
  
  names(rates.list)<-c(paste0(nsubset, "sub"))
  
  if(origin){
    cat("origin calculating\n")
    origin.aggr<-successTimes(collect, 500, 1, maxdim, mxscale, subsample = F)
    origin.rate<-lapply(1:maxdim, function(dim){
      return(successRates(origin.aggr, dim, corrects[dim]))
    })
    names(origin.rate)<-c(paste0("origine_dim", 1:maxdim))
    
    rates.list<-append(origin.rate, rates.list)
    
  }
  
  return(rates.list)
  
}


interpoCheck<-function(collect, nsubset, times, maxdim, mxscale, corrects, origin=F,
                       border=T){
  
  rates.list<-lapply(nsubset, function(nsub){
    cat(nsub, "subsample calculating\n")
    aggr.list<-interSuccessTimes(collect, nsub, times, maxdim, mxscale)
    rates<-lapply(1:maxdim, function(dim){
      rate<-successRates(aggr.list, dim, corrects[dim])
      inrate<-successRates(aggr.list, dim, corrects[dim], inter=T)
      return(append(rate, inrate))
    })
    cat(nsub, "subsample\n")
    names(rates)<-c(paste0("dim", 1:maxdim))
    print(rates)
    return(rates)
  })
  
  names(rates.list)<-c(paste0(nsubset, "sub"))
  
  if(origin){
    cat("origin calculating\n")
    origin.aggr<-interSuccessTimes(collect, 500, 1, maxdim, mxscale, subsample = F)
    origin.rate<-lapply(1:maxdim, function(dim){
      ori.rate<-successRates(origin.aggr, dim, corrects[dim])
      pri.inrate<-successRates(origin.aggr, dim, corrects[dim], inter=T)
      return(append(ori.rate, ori.inrate))
    })
    names(origin.rate)<-c(paste0("origine_dim", 1:maxdim))
    
    rates.list<-append(origin.rate, rates.list)
    
  }
  
  return(rates.list)
  
}

interSuccessTimes<-function(collect, nsub, times, maxdim, maxscale, subsample=T, border=T){
  
  aggr.list<-lapply(1:times, function(t){
    
    cat(t, "times subcollect\n")
    
    if(subsample){
      subcollect<-subsampleExclude(collect, nsub)
      aggr<-proposedMethodOnly.mk1(subcollect, maxdim, maxscale, 10)
      
      incollect<-subcollect
      for (l in 1:length(subcollect)) {
          inter.oricord<-voronoiInterpo(subcollect[[l]][["noizyX"]], 15)
          incollect[[l]][["noizyX"]]<-rbind(incollect[[l]][["noizyX"]], inter.oricord)
          incollect[[l]][["nsample"]]<-nrow(incollect[[l]][["noizyX"]])
          debugText(incollect[[l]][["nsample"]])
      }
      
      inaggr<-proposedMethodOnly.mk1(incollect, maxdim, maxscale, 10)
      
    }else{
      aggr<-proposedMethodOnly.mk1(collect, maxdim, maxscale, 10)
      incollect<-collect
      for (l in 1:length(collect)) {
        inter.oricord<-voronoiInterpo(subcollect[[l]][["noizyX"]], 15, border)
        incollect[[l]][["noizyX"]]<-rbind(incollect[[l]][["noizyX"]], inter.oricord)
        incollect[[l]][["nsample"]]<-nrow(incollect[[l]][["noizyX"]])
        debugText(incollect[[l]][["nsample"]])
      }
      
      inaggr<-proposedMethodOnly.mk1(subcollect, maxdim, maxscale, 10)
      
    }
    
    cycle<-lapply(1:maxdim, function(dim){
      
      result<-cyclenumber(aggr[[dim]], compare = F)
      cat(t, "times subcollect\n")
      cat(dim, "dimension\n")
      print(result)
      return(result)
      
    })
    
    incycle<-lapply(1:maxdim, function(dim){
      
      result<-cyclenumber(inaggr[[dim]], compare = F)
      cat(t, "times incollect\n")
      cat(dim, "dimension\n")
      print(result)
      return(result)
      
    })
    
    names(cycle)<-paste0("dim", 1:maxdim, "cycles")
    names(incycle)<-paste0("dim", 1:maxdim, "incycles")
    #debugText(incycle)
    return(c(aggr, cycle, inaggr, incycle))
    #return(append(aggr, cycle))
    
    
  })
  
  return(aggr.list)
  
}


####ボロノイ補間#######
voronoiProcess<-function(vics.line, figure){
  
  require(deldir)
  
  vics.pca<-prcomp(figure[vics.line,])
  
  res<-deldir(vics.pca$x[,1], vics.pca$x[,2])
  
  tiles<-tile.list(res)
  
  neibor1<-neighbourVoronoi(tiles, 1)
  
  neibor1<-c(1, neibor1)
  
  # for (i in neibor1) {
  #   
  #   ranpoint<-randomPointVoronoi(tiles[[i]])
  #   
  #   if(i==1){incord<-ranpoint}
  #   else{incord<-rbind(incord, ranpoint)}
  #   
  # }
  
  cenpoints<-sapply(neibor1, function(k)centerVoronoi(tiles[[k]]))
  
  vics.oricord<-originCoodinate(vics.pca, t(cenpoints))
  
  return(vics.oricord)
  
}

voronoiBorder<-function(vics.line, figure){
  
  require(deldir)
  
  vics.pca<-prcomp(figure[vics.line,])
  
  res<-deldir(vics.pca$x[,1], vics.pca$x[,2])
  
  tiles<-tile.list(res)
  
  insecs<-cbind(tiles[[1]][["x"]], tiles[[1]][["y"]])
  
  vics.oricord<-originCoodinate(vics.pca, insecs)
  
  return(vics.oricord)
  
}

voronoiInterpo<-function(figure, nvics, border=T){
  
  element<-rep(0, length = nrow(figure))
  
  dist<-distance(figure)
  
  for (i in 1:nrow(figure)) {
    
    if(element[i]==0){
      
      vics<-get.vicinity(dist, i, nvics)
      
      vics.line<-line.vics(i, vics)
      
      element[vics.line]<-element[vics.line]+1
      
      if(border==F){vics.oricord<-voronoiProcess(vics.line, figure)}
      else{vics.oricord<-voronoiBorder(vics.line, figure)}
      
      if(i==1){oricord<-vics.oricord}
      else{oricord<-rbind(oricord, vics.oricord)}
      
    }
    
  }
  
  #debugText(element)
  
  return(oricord)
  
}

centerVoronoi<-function(tile){
  
  cen.x<-mean(tile[["x"]])
  cen.y<-mean(tile[["y"]])
  
  return(c(cen.x, cen.y))
  
}
##########################################

#ネットから持ってきた関数
Print <- function(...){
  key <- as.list(substitute(list(...)))[-1L]
  val <- list(...)
  mapply(
    function(k, v){
      cat(k, "= ")
      if(!is.matrix(v) && (is.logical(v) || is.numeric(v) || is.complex(v) || is.character(v))){ cat(v, "\n") }
      else{ cat("\n"); print(v); cat("\n") }
    },  
    key, val)
  cat("\n")
}
#############################


#変数をRDataファイルに保存する関数
save2Rdata <- function(...) {
  require(tidyverse)
  elp <- list(...)
  elname <- substitute(...) %>% as.character()
  assign(elname, elp[[1]])
  save(list = elname, file = paste0("./data/", gsub("\\.", "_", elname), ".RData"))
}

#変数をRDataファイルに保存する関数2
#dataディレクトリがあるか確認し、なければ作成する
save2RData <- function(...) {
  require(tidyverse)
  elp <- list(...)
  elname <- substitute(...) %>% as.character()
  assign(elname, elp[[1]])
  
  dir<-match("data", list.files())
  
  if(is.na(dir)){
    
    dir.create("./data")
    
    
  }else{
    
    if(!(file.info("./data")$isdir)){dir.create("./data")}
  
    }
  
  save(list = elname, file = paste0("./data/", elname, ".RData"))
  
}
  

#補間を一括して行う関数
intering<-function(collect){
  incollect<-collect
  for (l in 1:length(collect)) {
    inter.oricord<-voronoiInterpo(collect[[l]][["noizyX"]], 15)
    incollect[[l]][["noizyX"]]<-rbind(incollect[[l]][["noizyX"]], inter.oricord)
    incollect[[l]][["nsample"]]<-nrow(incollect[[l]][["noizyX"]])
    debugText(l, incollect[[l]][["nsample"]])
  }
  
  return(incollect)
}

#補間の誤差計算
#トーラスに対して
#間違った関数
calc_error<-function(figure, maxr, minr, nps){
  
  require(tidyverse)
  
  error<-lapply((nps+1):(nrow(figure)), function(i){
    
    #debugText(i)
    
    al<-figure[i, 1]
    be<-figure[i, 2]
    
    x1<-sqrt(al^2+be^2) %>% {maxr*(al/.)}
    y1<-sqrt(al^2+be^2) %>% {maxr*(be/.)}
    
    dist<-dist(rbind(figure[i,], c(x1, y1, 0)))
    
    #debugText(dist)
    
    er<-abs(minr-dist)/minr
    
    #ebugText(er)
    
    return(er)
    
  })
  
  error<-unlist(error)
  
  return(error)
  
}

#補間の誤差計算
#トーラスに対して
torus_disterror<-function(figure, maxr, minr, nps){
  
  require(tidyverse)
  
  error<-lapply((nps+1):(nrow(figure)), function(i){
    
    #debugText(i)
    
    al<-figure[i, 1]
    be<-figure[i, 2]
    
    x1<-sqrt(al^2+be^2) %>% {maxr*(al/.)}
    y1<-sqrt(al^2+be^2) %>% {maxr*(be/.)}
    
    dist<-dist(rbind(figure[i,], c(x1, y1, 0)))
    
    #debugText(dist)
    
    dister<-dist-minr
    
    #ebugText(er)
    
    return(dister)
    
  })
  
  error<-unlist(error)
  
  return(error)
  
}

per_mean<-function(pd){
  
  mean<-zero_hat_double_threshold(pd)
  
  return(mean/2)
  
}


