#極値を探す関数
critical_point<-function(x){
  
  #値の差分をとる
  d<-diff(x)
  
  #極大値
  dec<-sapply(1:(length(d)-1), function(i){
    
    if(d[i]*d[i+1] <= 0 && d[i] >= 0){return(i+1)}
    
  }) %>% unlist()
  
  #極小値
  inc<-sapply(1:(length(d)-1), function(i){
    
    #debugText(i, x[i+1], x[i], d[i]*d[i+1])
    
    if(d[i]*d[i+1] <= 0 && d[i] < 0){return(i+1)}
    
  }) %>% unlist()
  
  dec<-dec[!is.null(dec)]
  inc<-inc[!is.null(inc)]
  
  #debugText(dec, inc)
  
  #隣接分を消す
  if(length(which(diff(dec)==1))>0){dec<-dec[-(which(diff(dec)==1)+1)]}
  if(length(which(diff(inc)==1))>0){inc<-inc[-(which(diff(inc)==1)+1)]}
  
  #端点
  if(d[1] < 0){dec<-c(1, dec)}
  if(d[length(d)] >= 0){dec<-c(dec, length(d)+1)}
  
  
  return(list(lmax=dec, lmin=inc))
  
}

#極値の組み合わせを求める
component_set<-function(crip){
  
  comp<-lapply(1:length(crip[["lmin"]]), function(i){
    
    return(c(crip[["lmax"]][i], crip[["lmin"]][i], crip[["lmax"]][i+1]))
    
  })
  
  return(t(as.data.frame(comp)))
  
}


#sublevel setのパーシステンスを調べる
sublevelset_persist<-function(x, crip){
  
  lmax<-order(x[crip$lmax]) %>% crip$lmax[.]
  lmin<-order(x[crip$lmin]) %>% crip$lmin[.]
  
  debugText(lmax, lmin)
  
  bdset<-matrix(0, length(crip$lmin), 2)
  
  for(i in length(crip$lmin):1){
    
    debugText(i)
    
    lmax2<-sort(lmax)
    
    debugText(lmin[i], lmax2)
    
    r_max<-lmax2[lmin[i] < lmax2][1]
    l_max<-lmax2[lmin[i] > lmax2][length(lmax2[lmin[i] > lmax2])]
    
    debugText(r_max, l_max)
    
    d<-which.min(x[c(r_max, l_max)]) %>% c(r_max, l_max)[.]
    
    debugText(d, lmax)
    
    if(d %in% lmax){
      
      bdset[i,]<-c(x[lmin[i]], x[d])
      
      lmax<-lmax[-which(lmax==d)]
      
      debugText(lmax, bdset)
      
    }
    
  }
  
  return(bdset)
  
}
