plot.result<-function(nset, dim, correct, Collect, Aggr){
  
  arrange<-matrix(0, nset, 7)
  for (p in 1:nset) {
    arrange[p, 1]<-Collect[[p]]$nsample
    arrange[p, 2]<-Collect[[p]]$var
    arrange[p, 3]<-Aggr[[dim]][p,"confidence"]
    arrange[p, 4]<-Aggr[[dim]][p,"bestconf"]
    arrange[p, 5]<-Aggr[[dim]][p,"landscape"]
    arrange[p, 6]<-Aggr[[dim]][p,"meanland"]
    arrange[p, 7]<-Aggr[[dim]][p,"proposed"]
  }
  
  colnames(arrange)<-c("nsample", "var", "confidence", "bestconf", "landscape", "meanland", "proposed")
  
  for(q in 1:5){
    plot(arrange[, 1], arrange[, 2],
         main = c("Confidence Interval", "Best Confidence Interval", "Persistent Landscape", "Mean Landscape", "Proposed Method")[q],
         xlab = "Data Points", ylab="Varience",
         col=ifelse((arrange[, q+2]>=(correct-0.5) & arrange[, q+2]<(correct+0.5)), "red", "blue"), type = "n")
    
    text(arrange[, 1], arrange[, 2],
         labels = round(arrange[, q+2], 1), col=ifelse((arrange[, q+2]>=(correct-0.5) & arrange[, q+2]<(correct+0.5)), "red", "blue"),
         cex = 0.8)
  }
  
  return(arrange)
  
}

plotWrongDiag<-function(Collect, Hole, method, correct){
  wrong<-which(Hole[, method]<(correct-0.5) | Hole[, method]>=(correct+0.5))
  
  #layout(matrix(1:((length(wrong)%/%4+1)*4), ncol=4))
  #split.screen(figs=c((length(wrong)%/%4)+1, 4))
 
  #pdf(paste0("./data/plot-wrong-diag-", method, ".pdf"))
  #dev.new()
  par(mfrow=c(((length(wrong)%/%4)+1), 4))
  for (k in 1:length(wrong)) {
    plot(Collect[[wrong[k]]]$diag)
  }
  #dev.off()
  par(mfrow=c(1, 1))
}

#分散を縦軸、データ点数を横軸に取り、推定サイクル数をプロット
#信頼区間、PL、MLと比較する場合とそれ以外で場合分け
plotAggr<-function(aggr, collect, correct, dim, compare=T, capture=0){
  nset<-length(collect)
  if(compare==T){
  arrange<-matrix(0, nset, 7)
    for (p in 1:nset) {
      arrange[p, 1]<-collect[[p]]$nsample
      arrange[p, 2]<-collect[[p]]$var
      arrange[p, 3]<-aggr[[dim]][p,"confidence"]
      arrange[p, 4]<-aggr[[dim]][p,"bestconf"]
      arrange[p, 5]<-aggr[[dim]][p,"landscape"]
      arrange[p, 6]<-aggr[[dim]][p,"meanland"]
      arrange[p, 7]<-aggr[[dim]][p,"proposed"]
    }
  
  colnames(arrange)<-c("nsample", "var", "confidence", "bestconf", "landscape", "meanland", "proposed")
  
  for(q in 1:5){
    plot(arrange[, 1], arrange[, 2],
         main = c("Confidence Interval", "Best Confidence Interval", "Persistent Landscape", "Mean Landscape", "Proposed Method")[q],
         xlab = "Data Points", ylab="Varience",
         col=ifelse((arrange[, q+2]>=(correct-0.5) & arrange[, q+2]<(correct+0.5)), "red", "blue"), type = "n")
    
    text(arrange[, 1], arrange[, 2],
         labels = round(arrange[, q+2], 1), col=ifelse((arrange[, q+2]>=(correct-0.5) & arrange[, q+2]<(correct+0.5)), "red", "blue"),
         cex = 0.8)
  }
  
  }else{
    arrange<-matrix(0, nset, 2+ncol(aggr[[dim]]))
    for (p in 1:nset) {
      arrange[p, 1]<-collect[[p]]$nsample
      arrange[p, 2]<-collect[[p]]$var
      for (q in 1:ncol(aggr[[dim]])) {
        arrange[p, q+2]<-aggr[[dim]][p,q]
      }
    }
    if(capture==0){colnames(arrange)<-c("nsample", "var", paste0(1:ncol(aggr[[dim]]), "method"))}
    else{colnames(arrange)<-c("nsample", "var", capture)}
    
    
    for(q in 1:ncol(aggr[[dim]])){
      plot(arrange[, 1], arrange[, 2],
           main = ifelse(capture==0, c(paste0("method",1:ncol(aggr[[dim]])))[q], capture[q]),
           xlab = "Data Points", ylab="Varience",
           col=ifelse((arrange[, q+2]>=(correct-0.5) & arrange[, q+2]<(correct+0.5)), "red", "blue"), type = "n")
      
      text(arrange[, 1], arrange[, 2],
           labels = round(arrange[, q+2], 1), col=ifelse((arrange[, q+2]>=(correct-0.5) & arrange[, q+2]<(correct+0.5)), "red", "blue"),
           cex = 0.8)
    }
    
  }
  
  return(arrange)
  
}

#リスト化された複数のPDを一括表示
plotPDs<-function(diags){
  
  par(mfrow=c(((length(diags)%/%4)+1), 4))
  for (k in 1:length(diags)) {
    plot(diags[[k]][[1]], cex.lab=1.6, cex.axis=1.6)
  }
  #dev.off()
  par(mfrow=c(1, 1))
}

#リスト化された複数のPDからPLを一括表示
plot_pls<-function(diags){
  
  par(mfrow=c(((length(diags)%/%4)*2+1), 4))
  lands<-lapply(diags, function(diag){
    
    land<-calcLandscape(diag)
    
  })
  #dev.off()
  par(mfrow=c(1, 1))
  
  return(lands)
  
}

#リスト化された2次PLを一括表示
plot_2ndpls<-function(lands, vert=F){
  
  oldpar<-par(no.readonly=T)
  #par(cex.lab=2, cex.main=2, cex.axis=2, plt = c(0.2, 0.9, 0.2, 0.9))
  
  if(vert==T){par(mfrow=c(2, (length(lands)%/%2+length(lands)%%2)), cex.lab=2, cex.main=2, cex.axis=2, plt = c(0.2, 0.9, 0.2, 0.9))}
  else{par(mfrow=c(((length(lands)%/%4)+1), 4), mgp=c(2.5, 1, 0))}
  for (k in 1:length(lands)) {
    plotLandscape(lands[[k]][-2])
  }
  #dev.off()
  #par(mfrow=c(1, 1))
  par(oldpar)

}
