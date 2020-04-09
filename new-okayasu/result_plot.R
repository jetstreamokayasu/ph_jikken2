library(MASS)

#平均、分散が等しいデータがパーシステントホモロジーで見分けられることを示す
mu0<-c(1, 1)

Sigma0<-rbind(c(1, 0.2), c(0.2, 1))

mv0<-mvrnorm(300, mu0, Sigma0)

##中心に穴
##400点生成して平均から近い順に200点減らす

mv1<-mvrnorm(400, mu0, Sigma0)

plot(mv1, xlab = "x", ylab = "y", xlim = c(-2.2, 4.2), ylim = c(-2.2, 4.2))

dist_mv1<-dist(rbind(mu0, mv1))

vanish_idx<-order(as.matrix(dist_mv1)[1,])[2:201]-1

vanished_mv1<-mv1[-vanish_idx, ]

mv2<-mvrnorm(400, mu0, Sigma0)

dist_mv2<-dist(rbind(mu0, mv2))

vanish_mv2_idx<-order(as.matrix(dist_mv2)[1,])[2:201]-1

vanished_mv2<-mv2[-vanish_mv2_idx, ]

plot(rbind(vanished_mv1, vanished_mv2), xlab = "x", ylab = "y", xlim = c(-2.2, 4.2), ylim = c(-2.2, 4.2))

x1<-rnorm(100)
y1<-rnorm(n = 100, mean = 0, sd = 1.5)

mv1_pd<-ripsDiag(mv1, maxdimension = 1, maxscale = 1.5)
vanished_pd<-ripsDiag(rbind(vanished_mv1, vanished_mv2), maxdimension = 1, maxscale = 1.5)