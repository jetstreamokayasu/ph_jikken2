source("new-okayasu/m1_main-mk1.R")
source("hikitsugi/Persistent Homology 計算とクイックルック.R")



# おかやす実験1 -----------------------------------------------------------------

nsample.min <- 200
nsample.max <- 600
var.min <- 0
var.max <- 2

ndata <- 2#データセット数
hole.set<-matrix(0, num.set, 4)
#サンプル数と分散、穴の推定数、穴の推定数の正しさを保存する配列
torus.500.test <- lapply(1:num.set,function(i){
  data.point<-round(runif(1, 200, 600))
  disper<-runif(1, 0,2)#分散
  hole.set[i,1]<<-data.point
  hole.set[i,2]<<-disper
  noize.torus <- matrix(0.3*rnorm(data.point*3, 0, disper), data.point)
  cat("i=", i, "data.point=", data.point, "disper=", disper, "\n")
  return(torusUnif(data.point, 1,2.5) + noize.torus)
})

data.collect <- lapply(1:ndata, function(i){
  nsample <- round(runif(1, nsample.min, nsample.max))
  var <- runif(1, var.min, var.max)
  noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, var = var, noizyX = torus + noize.torus,
              X = torus, diag = 0))
})

for (i in 1:ndata) {
  data.collect[[i]]$diag <- 123
}



set.seed(1908)
ex1.test.aggr <- homologyMethodsComp(torus.500.test,2,3,hole.set[,1],10)
#round(data.point*(3/5))
#サンプル数をランダムに設定する場合、
#データセットごとにサブサンプル数を変更する必要がある？

#暫定処置としてサブサンプル数を300に固定し、
#サンプル数は300~600の一様分布からのランダムとする

Discrimination<-function(X){
  for (j in 1:num.set) {
    
    if(X[j, 3]==2) X[j, 4]<-2
    #穴の推定数が正しければ2を代入
    
    else X[j, 4]<-4#穴の推定数が正くなければければ4を代入
    
  }
  
  plot(X[,1], X[,2], col=X[,4])
}

Discrimination(hole.set)
