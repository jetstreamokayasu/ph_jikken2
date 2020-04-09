torus.500.10 <- lapply(1:10,function(i)torusUnif(500,1,2.5))
set.seed(1908)
ex1.aggr <- homologyMethodsComp(torus.500.10,2,3,300,10)



torus.500.100 <- lapply(1:100,function(i)torusUnif(500,1,2.5))
for(i in 1:10){
  set.seed(1908)
  ex2.temp <- homologyMethodsComp(torus.500.100[((i-1)*10+1):(i*10)],2,3,300,10)
  if(i==1) ex2.aggr <- ex2.temp
  else ex2.aggr <- lapply(1:2,function(d)ex2.aggr[[d]]+ex2.temp[[d]])
}
# set.seed(1908)
# ex2.aggr <- homologyMethodsComp(torus.500.100,2,3,300,10)

# set.seed(1908)
# ex3.aggr <- homologyMethodsComp(torus.500.100,2,3,size = 300,samples = 30)

sphere.300.10 <- lapply(1:10,function(i)sphereUnif(300,d = 2))
set.seed(1908)
ex3.aggr <- homologyMethodsComp(sphere.300.10,2,1.9,size = 200,samples = 10)

sphere.500.100 <- lapply(1:100,function(i)sphereUnif(500,2))
set.seed(1908)
ex4.aggr <- homologyMethodsComp(sphere.500.100,2,1.9,size = 300,samples = 10)


# Centroid Thresh ---------------------------------------------------------
set.seed(1908)
ex5.aggr <- homologyMethodsComp(torus.500.100,2,3,300,10)

set.seed(1908)
ex6.aggr <- homologyMethodsComp(torus.500.100,2,3,300,10)


# Centroid times sqrt(2) /2 thresh ----------------------------------------
set.seed(1908)
ex7.aggr <- homologyMethodsComp(torus.500.100,2,3,300,10)


sphere.300.100 <- lapply(1:100,function(i)sphereUnif(300,2))
set.seed(1908)
ex7.aggr <- homologyMethodsComp(sphere.300.100,2,1.9,size = 200,samples = 10)

hausdInterval(torus.500.100[[1]],499,30)
