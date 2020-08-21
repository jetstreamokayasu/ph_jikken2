#高次元データを試す
#3次元球面から
#s(n)phereのnの数で次元を表す
s4phere<-sphereUnif(100, 3)
s3phere<-sphereUnif(100, 2)

s4phere_pd<-ripsDiag(X = s4phere, maxdimension = 4, maxscale = 2, printProgress = T)
