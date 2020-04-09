 #start, endで区切られた範囲にデータ点が存在しするかしないかを判定する関数
#存在する場合はT、存在しない場合はFを返す
existCheck<-function(start, end, mapped){
  
  exie1.mem<-which(mapped[,1]>=start[1] & mapped[,1]<end[1])
  
  #debugText(exie1.mem, start, end)
  
  if(length(exie1.mem)>=1){
    
    exie2.mem<-which(mapped[exie1.mem, 2]>=start[2] & mapped[exie1.mem ,2]<end[2])
    
    #debugText(exie2.mem)
    
    if(length(exie2.mem)>=1) return(T)
    
    else return(F)
    
  }
  
  else return(F)
  
}

#格子状の線を引く
gridLine<-function(x, div){
  
  xlim<-range(x[,1])
  ylim<-range(x[,2])
  
  wid<-abs(xlim[2]-xlim[1])/div
  hei<-abs(ylim[2]-ylim[1])/div
  
  for (i in 0:div) {
    
    abline(v=xlim[1]+i*wid)
    abline(h=ylim[1]+i*hei)
    
  }
  
}


pixelConvert<-function(x, div){
  
  xlim<-range(x[,1])
  ylim<-range(x[,2])
  
  wid<-abs(xlim[2]-xlim[1])/div
  hei<-abs(ylim[2]-ylim[1])/div
  
  pixel<-matrix(0, div, div)
  
  for (i in 1:div) {
    
    for (j in 1:div) {
      
      if(existCheck(c(xlim[1]+(i-1)*wid, ylim[1]+(j-1)*hei), c(xlim[1]+i*wid, ylim[1]+j*hei), x)){
    
            pixel[j, i]<-1
        
      }
    }
    
  }
  
  return(pixel)
    
  }

#行列の注目した要素の8近傍に要素が入っているか
#入っていればTを返す
neighEleCheck<-function(pic, row, col){
  
  if(row!=1 && (col!=1) && (pic[row-1, col-1]>0)){return(T)}
  else if((row!=1) && pic[row-1, col]>0){return(T)}
  else if((row!=1) && col!=ncol(pic) && pic[row-1, col+1]>0){return(T)}
  else if(col!=1 && pic[row, col-1]>0){return(T)}
  else if(col!=ncol(pic) && pic[row, col+1]>0){return(T)}
  else if(row!=nrow(pic) && col!=1 && pic[row+1, col-1]>0){return(T)}
  else if(row!=nrow(pic) && pic[row+1, col]>0){return(T)}
  else if(row!=nrow(pic) && col!=ncol(pic) && pic[row+1, col+1]>0){return(T)}
  else{return(F)}
  
}  

#注目要素に要素が無く、注目要素の8近傍に要素があるとき
#注目要素に2を代入
insertElement<-function(pic){
  
  cp.pic<-pic
  
  for (i in 1:nrow(pic)) {
    
    for (j in 1:ncol(pic)) {
      
      if(pic[i,j]==0 &&neighEleCheck(pic, i, j)){
        
        cp.pic[i, j]<-2
        
      }
      
    }
    
  }
  
  return(cp.pic)
  
}

#指定された範囲にPCA後の座標を返す
pcaCoordinate<-function(xmin, ymin, wid, hei, row, col){
  
  return(c(xmin+(wid*(col-1))+wid/2, ymin+(hei*(row-1))+hei/2))
  
}

#データ点が存在しないピクセルの中央にデータ点を打つ
#座標はPCA後の座標
pcaCoord.set<-function(x, cppic, div){
  
  xlim<-range(x[,1])
  ylim<-range(x[,2])
  
  wid<-abs(xlim[2]-xlim[1])/div
  hei<-abs(ylim[2]-ylim[1])/div
  
  ele<-which(cppic==2, arr.ind=TRUE)
  
  coord<-sapply(1:nrow(ele), function(k)pcaCoordinate(xlim[1], ylim[1], wid, hei, ele[k,1], ele[k,2]))
  
  return(t(coord))
  
}

#PCA適用後の座標から元座標を算出
originCoodinate<-function(rpca, incord){
  
  eigen01<-as.matrix(rpca$rotation[,1])
  eigen02<-as.matrix(rpca$rotation[,2])
  
  if(!is.matrix(incord)){incord<-t(as.matrix(incord))}
  
  oricord<-sapply(1:nrow(incord), function(l){
    
  return((incord[l, 1]*(eigen01)+incord[l, 2]*(eigen02))+rpca$center)
    
  })
  
  return(t(oricord))
  
}

#任意の点に近傍点の中から最も遠い点を抽出
#そのもっとも遠い点から近傍の中で最も近いn点を抽出
coveredVic<-function(vicsline, figure, n){
  
  distfarvic<-sapply(1:(length(vicsline)-2), function(k){
    
    dist.set<-c(vicsline[length(vicsline)], vicsline[k+1], sum((figure[vicsline[k+1],]-figure[vicsline[length(vicsline)],])^2))
    names(dist.set)<-c("start", "goal", "distance")
    
    return(dist.set)
    
  })
  
  distfarvic<-t(distfarvic)
  
  distfarvic<-distfarvic[order(distfarvic[,3]),]
  
  return(distfarvic[1:n,])
  
}

#指定されたデータ点のnvics点近傍をPCAで変換し
#膨張処理を行い、補間されたデータ点の下\元の座標系での座標を返す
expandProcess<-function(vics, vics.line, figure, dist, div){
  
  vics.pca<-prcomp(figure[vics.line,])
  
  vics.pic<-pixelConvert(vics.pca[["x"]], div)
  
  vics.cppic<-insertElement(vics.pic)
  
  vics.incord<-pcaCoord.set(vics.pca[["x"]], vics.cppic, div)
  
  vics.oricord<-originCoodinate(vics.pca, vics.incord)
  
  return(vics.oricord)
  
}

interPolation_test<-function(figure, nvics, div){
  
  element<-rep(0, length = nrow(figure))
  
  dist<-distance(figure)
  
  for (i in 1:nrow(figure)) {
    
    if(element[i]==0){
      
      vics<-get.vicinity(dist, i, nvics)
      
      vics.line<-line.vics(i, vics)
      
      element[vics.line]<-element[vics.line]+1
      
      vics.oricord<-expandProcess(vics, vics.line, figure, dist, div)
      
      if(i==1){oricord<-vics.oricord}
      else{oricord<-rbind(oricord, vics.oricord)}
      
    }
    
  }
  
  #debugText(element)
  
  return(oricord)
  
}

#3次元図形をプロット
figurePlot<-function(X){
  
  require(rgl)
  plot3d(X)
  aspect3d("iso")
  
}

distance<-function(origin){
  
  dist<-matrix(0, (1/2)*nrow(origin)*(nrow(origin)-1), 3)
  colnames(dist)<-c("start", "goal", "distance")
  r<-1
  
  for (k in 1:nrow(origin)) {
    for (l in (k+1):nrow(origin)) {
      
      #if((k!=l) && !(k %in% dist[,1]) && !(l %in% dist[,2])){
      if(l <= nrow(origin) && k!=l){
        #debugText(k ,l, r)
        dist[r,1]<-k
        dist[r,2]<-l
        dist[r,3]<-sum((origin[k,]-origin[l,])^2)
        r<-r+1
      }
      #}
    }
    
  }
  return(dist)
}

#任意(center)の点から最も近いnvic点を求める関数
get.vicinity<-function(dis, center, nvic){
  
  choice<-rbind(dis[which(dis[, "start"]==center), ], dis[which(dis[, "goal"]==center), ])
  choice<-choice[order(choice[,3]),]
  
  vic<-choice[1:nvic,]
  
  return(vic)
  
}

#ある点を中心とした近傍点を色付けてプロット
figurePlot.coloredVic<-function(figure, vics, centr){
  
  require(rgl)
  
  vics.line<-line.vics(centr, vics)
  
  plot3d(figure[-vics.line, ])
  aspect3d("iso")
  points3d(figure[vics.line, ], col=3)
  
}

line.vics<-function(centr, vic){
  
  vics<-sapply(1:length(vic[,2]), function(t){
    
    if(vic[t,2]==centr) return(vic[t, 1])
    else return(vic[t, 2])
    
  })
  
  names(vics)<-NULL
  
  vics<-c(centr, vics)
  
  return(vics)
  
}

conbineInterOrigin<-function(figure, interpo){
  
  return(rbind(figure, interpo))
  
}

meanVicsDestance<-function(figure, nvics){
  
  element<-rep(0, length = nrow(figure))
  
  dist<-distance(figure)
  
  for (i in 1:nrow(figure)) {
    
    if(element[i]==0){
      
      vics<-get.vicinity(dist, i, nvics)
      
      if(i==1){vics.maxdist<-vics[nvics, "distance"]}
      else{vics.maxdist<-c(vics.maxdist, vics[nvics, "distance"])}
      
      vics.line<-line.vics(i, vics)
      
      element[vics.line]<-element[vics.line]+1
      
    }
    
  }
  
  #debugText(element)
  
  return(mean(vics.maxdist))
  
}

meanInterPolation<-function(figure, nvics){
  
  element<-rep(0, length = nrow(figure))
  
  dist<-distance(figure)
  
  for (i in 1:nrow(figure)) {
    
    if(element[i]==0){
      
      vics<-get.vicinity(dist, i, nvics)
      
      vics.line<-line.vics(i, vics)
      
      element[vics.line]<-element[vics.line]+1
      
      mean.vics<-apply(figure[vics.line,], 2, mean)
      
      if(i==1){inter<-mean.vics}
      else{inter<-rbind(inter, mean.vics)}
      
    }
    
  }
  
  #debugText(element)
  
  return(inter)
  
}

#ボロノイ領域のある頂点を含む辺を出力
vertex.side<-function(tile, vertex){
  
  vert.set<-1:length(tile[["x"]])
  vert.set<-c(length(tile[["x"]]), vert.set, 1)
  
  sides<-matrix(0, 2, 2)
  sides[1,]<-vert.set[c(vertex, vertex+1)]
  sides[2,]<-vert.set[c(vertex+1, vertex+2)]
  
  return(sides)
  
}

#凸包内のある頂点を含む辺を出力
convex_hull_vertx<-function(chul, vertx){
  
  
  vert.set<-c(chul[length(chul)], chul, chul[1])
  
  #debugText(vert.set, which(chul==vertx))
  
  sides<-matrix(0, 2, 2)
  sides[1,]<-which(chul==vertx) %>% '+'(., c(0, 1)) %>%   vert.set[.]
  sides[2,]<-which(chul==vertx) %>% '+'(., c(1, 2)) %>%   vert.set[.]
  
  return(sides)
  
}

#PCAで写された点が作る凸包内に補間点があるか判定
convex_hull_check<-function(rpca, hline, sides){
  
  t1<-sapply(sides[,1], function(side){
    
    return((hline[1,1]-hline[2,1])*(rpca[["x"]][side, 2]-hline[1,2]))
    
  })
  
  t2<-sapply(sides[,2], function(side){
    
    return((hline[1,1]-hline[2,1])*(rpca[["x"]][side, 2]-hline[1,2]))
    
  })
  
  
  t3<-sapply(1:nrow(sides), function(k){
    
    return((rpca[["x"]][sides[k,1], 1]-rpca[["x"]][sides[k,2], 1])*(hline[1,2]-rpca[["x"]][sides[k,1], 2])+(rpca[["x"]][sides[k,1], 2]-rpca[["x"]][sides[k,2], 2])*(rpca[["x"]][sides[k,1], 1]-hline[1,1]))
    
  })
  
  
  t4<-sapply(1:nrow(sides), function(k){
    
    return((rpca[["x"]][sides[k,1], 1]-rpca[["x"]][sides[k,2], 1])*(hline[2,2]-rpca[["x"]][sides[k,1], 2])+(rpca[["x"]][sides[k,1], 2]-rpca[["x"]][sides[k,2], 2])*(rpca[["x"]][sides[k,1], 1]-hline[2,1]))
    
  })
  
  #debugText(t1, t2, t3, t4)
  # debugText((t1*t2)<0)
  # debugText((t3*t4)<0)
  #ncross<-length(which((t1*t2)<0))
  ncross<-((t1*t2)<0 & (t3*t4)<0)
  
  return(ncross)
  
}

#渡された点は凸包内にあるか判定
exist_convexhull_check<-function(rpca, insecs){
  
  chul<-chull(rpca[["x"]][,1:2])
  
  exist<-sapply(1:nrow(insecs), function(i){
    sides<-chul[which(rpca[["x"]][chul,1]>=insecs[i,1])] %>% 
                sapply(., function(k)convex_hull_vertx(chul, k))
    if(length(sides)==0){return(F)}
    else{cross.side<-sidesSet(sides)}
    hline<-matrix(c(insecs[i,], max(rpca[["x"]][chul,1][which(rpca[["x"]][chul,1]>=insecs[i,1])]), insecs[i,2]), 2, 2, byrow=T)
    #debugText(hline)
    c.ncross<-convex_hull_check(rpca, hline, t(cross.side))
    #debugText(c.ncross)
    if(length(which(c.ncross==T)) %% 2 != 0){return(T)}
    else{return(F)}
    
  })
  
  return(exist)
  
}

#交差判定を行う辺の頂点集合を作成
# vertexSet<-function(tile, sides){
#   
#   if(ncol(sides)==1){return(sort(sides[,1]))}
#   
#   else{
#   for (i in 1:(ncol(sides)-1)) {
#     debugText(i)
#     if(i==1){ver.set<-union(sides[,i], sides[,i+1])}
#     else{ver.set<-union(ver.set, sides[,i+1])}
#     debugText(ver.set)
#   }
#   
#   return(sort(ver.set))
#   }
#   
# }

#辺の集合で被りを無くす
sidesSet<-function(sides){
  #debugText(sides)
  if(!is.matrix(sides)){sides<-t(as.matrix(sides))}
  #debugText(sides)
  check.sides<-matrix(0, 2, ncol(sides)*2)
  #debugText(ncol(sides))
  t<-1
  
  for (i in 1:ncol(sides)) {
    check.sides[,t]<-c(sides[1,i], sides[2,i])
    t<-t+1
    check.sides[,t]<-c(sides[3,i], sides[4,i])
    t<-t+1
  }
  
  for(j in seq(ncol(check.sides), 1)){
    
    #debugText(j)
    
    for(k in seq((ncol(check.sides)-1), 1)){
      
      if(j!=k && setequal(check.sides[,j], check.sides[,k])){
        
        #debugText(k)
        
        check.sides<-check.sides[,-j]
        
        #debugText(check.sides)
        
        break
        
      }
      
    }
    
  }
  
  return(check.sides)
  
}

#ボロノイ領域内に点があるか調べる
crossCheck<-function(tile, hline, sides){
  
  t1<-sapply(sides[,1], function(side){
    #debugText(side)
    
    return((hline[1,1]-hline[2,1])*(tile[["y"]][side]-hline[1,2]))
    
  })
  
  t2<-sapply(sides[,2], function(side){
    
    #cat("t2_side=", side, "\n")
    
    return((hline[1,1]-hline[2,1])*(tile[["y"]][side]-hline[1,2]))
    
  })

  
  t3<-sapply(1:nrow(sides), function(k){
    
    return((tile[["x"]][sides[k,1]]-tile[["x"]][sides[k,2]])*(hline[1,2]-tile[["y"]][sides[k,1]])+(tile[["y"]][sides[k,1]]-tile[["y"]][sides[k,2]])*(tile[["x"]][sides[k,1]]-hline[1,1]))
    
  })
  
  
  t4<-sapply(1:nrow(sides), function(k){
    
    return((tile[["x"]][sides[k,1]]-tile[["x"]][sides[k,2]])*(hline[2,2]-tile[["y"]][sides[k,1]])+(tile[["y"]][sides[k,1]]-tile[["y"]][sides[k,2]])*(tile[["x"]][sides[k,1]]-hline[2,1]))
    
  })
  
  # debugText(t1, t2, t3, t4)
  # debugText((t1*t2)<0)
  # debugText((t3*t4)<0)
  
  #ncross<-length(which((t1*t2)<0))
  ncross<-((t1*t2)<0 & (t3*t4)<0)
  
  return(ncross)
  
}

#あるボロノイ領域内にランダムに点を打つ
randomPointVoronoi<-function(tile){
  
  inter<-F
  while (inter==F) {
    
    ranx<-runif(1, range(tile[["x"]]))
    rany<-runif(1, range(tile[["y"]]))
  
    cross.mem<-which(tile[["x"]]>=ranx)
    sides<-sapply(cross.mem, function(k)vertex.side(tile, k))
    if(length(sides)<1){next}
    check.side<-sidesSet(sides)
    hline<-matrix(c(ranx, rany, max(tile[["x"]][which(tile[["x"]]>=ranx)]), rany), 2, 2, byrow=T)
    ncross<-crossCheck(tile, hline, t(check.side))
    
    if(length(which(ncross==T)) %% 2 != 0){inter=T}
    
  }
  
  return(c(ranx, rany))
  
}

#隣接するボロノイ領域を探す
neighbourVoronoi<-function(tiles, centr){
  
  neibor<-sapply(1:length(tiles), function(k){
    
    if(k!=centr && any(tiles[[centr]][["x"]] %in% tiles[[k]][["x"]]))
      
      return(k)
    
  })
  
  for (i in length(neibor)) {
    
    if(is.null(neibor[[i]])){neibor[i]<-NULL}
    
  }
  
  return(unlist(neibor))
    
}

voronoiInterpo<-function(figure, nvics){
  
  element<-rep(0, length = nrow(figure))
  
  dist<-distance(figure)
  
  for (i in 1:nrow(figure)) {
    
    if(element[i]==0){
      
      vics<-get.vicinity(dist, i, nvics)
      
      vics.line<-line.vics(i, vics)
      
      element[vics.line]<-element[vics.line]+1
      
      #vics.oricord<-voronoiProcess(vics.line, figure)
      vics.oricord<-voronoiBorder(vics.line, figure)[[1]]
      
      if(i==1){oricord<-vics.oricord}
      else{oricord<-rbind(oricord, vics.oricord)}
      
    }
    
  }
  
  #debugText(element)
  
  return(oricord)
  
}

#指定されたデータ点のnvics点近傍をPCAで変換し
#ボロノイ図を描き、中心点のボロノイ領域
#及び隣接するボロノイ領域内にランダムに点を打ち補間
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

#ボロノイ領域の中心に点を打つ
centerVoronoi<-function(tile){
  
  cen.x<-mean(tile[["x"]])
  cen.y<-mean(tile[["y"]])
  
  return(c(cen.x, cen.y))
  
}

#ボロノイ領域の頂点に点を打つ
#PCAで写された点による凸包内に入っていない点は除く
voronoiBorder<-function(vics.line, figure){
  
  require(deldir)
  
  vics.pca<-prcomp(figure[vics.line,])
  
  res<-deldir(vics.pca$x[,1], vics.pca$x[,2])
  
  tiles<-tile.list(res)

  insecs<-cbind(tiles[[1]][["x"]], tiles[[1]][["y"]])
  
  exist<-exist_convexhull_check(vics.pca, insecs)
  
  #debugText(vics.line, exist, insecs[which(exist==T), ])
  
  vics.oricord<-originCoodinate(vics.pca, insecs[which(exist==T), ])
  
  return(list(oricord=vics.oricord, pca.inter=insecs[which(exist==T), ]))
  
}

#補間した点の誤差を求める
#求まっていない。要修正
errorTorus<-function(r, R, ori, intered){
  
 errors<-sapply((ori+1):length(intered[,1]), function(k){
   if(abs(intered[k,3])>r){
     
     return(abs(r-abs(intered[k,3]))/r)
     #cat(k, "=", 1, "\n")
     }
   else{
     
     sqxyR<-(sqrt(intered[k,1]^2+intered[k,2]^2)-R)^2
     
     if(sqxyR>r^2){
       
       error<-abs(r-sqxyR)/r 
       #cat(k, "=", 2, "\n")
       return(error)
       
       }
     
     else{
       
       truth.z<-sqrt(r^2-sqxyR)
       #print(sqxyR)
       error<-(abs(truth.z-abs(intered[k,3])))/(truth.z)
       
       #debugText(intered[k,3], truth.z, error)
       #cat(k, "=", 3, "\n")
       
       return(error)
       
     }
     
   }
   
 })
 
 return(errors*100)
  
}

#PCAで求まった基底によって張られた平面の方程式の係数を求める
#centrは元の多様体(figure)の近傍ｋ点の中心点。元の多様体の何番目の点かで指定
confirmPlane<-function(rpca){
  
  cros<-pracma::cross(rpca[["rotation"]][,1], rpca[["rotation"]][,2])
  
  seg<-originCoodinate(rpca, t(as.matrix(rpca[["x"]][1, 1:2])))
  
  coefs<-c(cros, -sum(cros*seg))
  
  return(coefs)
  
}
