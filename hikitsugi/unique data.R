library(TDA)
library(rgl)

appleUnif <- function(n,depth=2.5){
  apple <- sphereUnif(n,2)
  apple <- apple * sin(apple[,3]*depth)
  apple
}

swissUniff <- function(n,curveness=2,length=2,width=4){
  if(curveness<length) stop("curveness must be equal or larger than length")
  p = sqrt(curveness + length * seq(-1, 1 - 2/n, 2/n))
  y = width/2 * runif(n, -1, 1)
  d_sr = cbind(p * cos(2*pi*p), y, p * sin(2*pi*p))
  d_sr
}

anulusUnif <- function(n,r.in=1,r.out=2){
  theta <- runif(n,0,2*pi)
  r <- sqrt(2*runif(n)/(2/(r.out^2-r.in^2))+r.in^2)
  cbind(r*cos(theta),r*sin(theta))
}

vonMisesUnif <- function(){
  
}

cylinderUnif <- function(n,height=1){
  theta <- runif(n,min = 0,max = 2*pi)
  l <- runif(n,min = -height/2,height/2)
  cbind(cos(theta),sin(theta),l)
}

scurveUnif <- function(n,width=1){
  t <- 3*pi*runif(n,-0.5,0.5)
  x <- sin(t)
  y <- width*runif(n)
  z <- sign(t) * (cos(t)-1)
  cbind(x,y,z)
}

sanulusUniff <- function(n,r.in=1,r.out=2){
  an <- anulusUnif(n,r.in = r.in,r.out = r.out)
  f <- function(x) ifelse(x>0,-3/4*pi*x+2*pi,3/4*pi*x+pi)
  t <- f(an[,1])*sign(an[,1])
  x <- sin(t)*sign(an[,1])
  z <- (cos(t)-1) - sign(an[,1])
  cbind(x,an[,2],z)
}

xplaneUnif <- function(n,r.in,r=1,r.out=2){
  t <- 4*pi*runif(n,-0.5,0.5)
  x <- cos(t)
  y <- sin(t)
  z <- sign(t) * (cos(t*2)-1)
  cbind(x,y,z)
}

xanulusUnif <- function(n,r.in=1,r.out=2){
  t <- 4*pi*runif(n,-0.5,0.5)
  r <- sqrt(2*runif(n)/(2/(r.out^2-r.in^2))+r.in^2)
  x <- cos(t)*r
  y <- sin(t)*r
  z <- sign(t) * (cos(t*2)-1)
  cbind(x,y,z)
}

torsionScurveUnif <- function(n,r.in=1,r.out=2){
  t <- 3*pi*runif(n,-0.5,0.5)
  r <- sqrt(2*runif(n)/(2/(r.out^2-r.in^2))+r.in^2)
  x <- sign(t)*r-sign(t)*1.5
  y <- sin(t)*r
  z <- sign(t) * (cos(t)-1)
  cbind(x,y,z)
}

strangeAnulusUnif <- function(n,r.in=1,r.out=2){
  t <- 2*pi*runif(n,-0.5,0.5)
  r <- sqrt(2*runif(n)/(2/(r.out^2-r.in^2))+r.in^2)
  x <- cos(t)*r
  y <- sin(t)*r
  z <- sin(y) * (cos(t)) * sign(t)
  cbind(x,y,z)
}

torsionAnulusUnif <- function(n,r.in=1,r.out=2){
  t <- 2*pi*runif(n,-0.5,0.5)
  r <- sqrt(2*runif(n)/(2/(r.out^2-r.in^2))+r.in^2)
  x <- cos(t)*r
  y <- sin(t)*r
  z <- sin(y) * (cos(t))
  cbind(x,y,z)
}

twistAnulusUnif <- function(n,r.in=1,r.out=2){
  t <- 2*pi*runif(n,-0.5,0.5)
  r <- sqrt(2*runif(n)/(2/(r.out^2-r.in^2))+r.in^2)
  x <- cos(t)*r * sin(t)^2
  y <- sin(t)*r
  z <- sin(y) * (cos(t))
  cbind(x,y,z)
}

aHoleRollUnif <- function(n,r.in=1,r.out=2){
  an <- anulusUnif(n,r.in,r.out)
  f <- function(x) 3/4*pi*x + pi/2
  t <- f(an[,1])*sign(an[,1])
  x <- sin(t)
  z <- sign(t) * (cos(t)-1)
  cbind(x,an[,2],z)
}

s3dUnif <- function(n,r.in=1,r.out=2){
  theta <- runif(n,0,2*pi)
  t <- 3*pi*runif(n,-0.5,0.5)
  r <- sqrt(2*runif(n)/(2/(r.out^2-r.in^2))+r.in^2)
  x <- r*cos(t)*sin(t)
  y <- r*sin(t)
  z <- sign(t)*(cos(t)-1)
  cbind(x,y,z)
}

holledRollUnif <- function(n,r.in=1,r.out=2){
  an <- anulusUnif(n,r.in,r.out)
  f <- function(x) 3/4*pi*x + pi/2
  t <- f(an[,1])*sign(an[,1])+(an[,1]>0*pi)
  hist(t/pi)
  x <- sin(t)
  z <- sign(t) * (cos(t)-2)
  cbind(x,an[,2],z)
}

helixAnulusUnif <- function(n,spin = 1){
  an <- anulusUnif(n)
  theta <- an[,1] * pi / 2 * spin
  return(cbind(an[,1],cos(theta)*an[,2],sin(theta)*an[,2]))
}

anulusUnif <- function(n,r.in=1,r.out=2){
  theta <- sort(runif(n,0,2*pi))
  r <- sqrt(2*runif(n)/(2/(r.out^2-r.in^2))+r.in^2)
  cbind(r*cos(theta),r*sin(theta))
}

# spiral ------------------------------------------------------------------
cylinder <- function(theta,l){
  x <- cos(theta)
  y <- sin(theta)
  z <- l
  cbind(x,y,z)
}

anulusReg <- function(n,r.in=1,r.out=2){
  theta <- seq(0,2*pi,length.out = n)
  r <- sqrt(2*runif(n)/(2/(r.out^2-r.in^2))+r.in^2)
  cbind(r*cos(theta),r*sin(theta)) 
}

panelReg <- function(n,length=2,width=1){
  x <- runif(n,min = -length/2,max = length/2)
  y <- runif(n,min = -width/2,max = width/2)
  return(cbind(sort(x),y))
}

#                1,function(x)cylinder(x[1],x[2]))),col=rainbow(20)[seq(0,1,length.out = 20)*30])