#アクロボット画像もどきを作成
theta<-seq(0, 2*pi, length=20)-pi/2
plot(cos(theta), sin(theta), xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), type="n", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
text(cos(theta), sin(theta), 1:20)
lines(c(0, cos(theta[6])), c(0, sin(theta[6])), lwd=10)

alpha<-seq(0, 2*pi, length=20)
points(0.5*cos(alpha)+cos(theta[6]), 0.5*sin(alpha)+sin(theta[6]))
lines(c(cos(theta[6]), 0.5*cos(alpha[6])+cos(theta[6])), c(sin(theta[6]), 0.5*sin(alpha[6])+sin(theta[6])), lwd=10)

oldpar <- par(no.readonly = TRUE) 
par(mar=c(0.1, 0.1, 0.1, 0.1), omi = c(0.1, 0.1, 0.1, 0.1), xpd=NA)
par(oldpar)


{
k<-1
r1<-1
r2<-0.5
a<-0.1

for(i in 1:length(theta)){
  
  for(j in 1:length(alpha)){
    
    png(filename=paste0("./data/acro_img2/", k, ".png"), width = 50, height = 50)
    #par(mar=c(0, 0, 0, 0), omi = c(0, 0, 0, 0))
    #par(plt = c(0, 1, 0, 1), omd = c(0, 1, 0, 1))
    par(mar=c(0.1, 0.1, 0.1, 0.1), omi = c(0.1, 0.1, 0.1, 0.1))
    #plot(r1*cos(theta), r1*sin(theta), xlim = c(-(r1+r2+a), (r1+r2+a)), ylim = c(-(r1+r2+a),(r1+r2+a)))
    plot(r1*cos(theta), r1*sin(theta), xlim = c(-(r1+r2+a), (r1+r2+a)), ylim = c(-(r1+r2+a), (r1+r2+a)), type="n", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
    lines(c(0, r1*cos(theta[i])), c(0, r1*sin(theta[i])), lwd=5)
    lines(c(r1*cos(theta[i]), r2*cos(alpha[j])+r1*cos(theta[i])), c(r1*sin(theta[i]), r2*sin(alpha[j])+r1*sin(theta[i])), lwd=5)
    dev.off()
    k<-k+1
  }
}
}

files<-list.files("./data/acro_img2/")
file.remove("\\.png$", paste0("./data/acro_img2/", files))


#アクロボット画像の分析
require(png)
require(myimg)
require(rgl)
arrow <- lapply(paste0(1:(arrow.row * arrow.col), ".png"), readPNG)
arrow <- lapply(arrow, convert2Gray)
arrow.mat <- matrix(unlist(arrow), arrow.row * arrow.col, 50 * 25, byrow = T)
arrow.dist <- dist(arrow.mat)

acro2<-lapply(paste0("./data/acro_img2/",1:400, ".png"), readPNG)
acro2_grays<-lapply(acro2, convert2Gray)
showImage(acro2_grays[[10]])
showImageList(acro2_grays[seq(20, 100, 10)+0:8])
acro2_mat<-matrix(unlist(acro2_grays), 400, 50*50, byrow = T)

acro2_pca<-prcomp(acro2_mat)
plot3d(acro2_pca[["x"]][,1:3])

acro2_dist<-dist(acro2_mat)
acro2_diag<-ripsDiag(acro2_mat, maxdimension = 2, maxscale = 10)
plot(acro2_diag[[1]])

# lle
acro_lle <- lle::lle(acro_mat, m = 3,k = 12)$Y
plot3d(acro_lle, col = rainbow(400))

# isomap
arrow.iso <- vegan::isomap(arrow.dist, ndim = 3, k=12)[[1]]
plot3d(arrow.iso, col = rainbow(arrow.row * arrow.col))