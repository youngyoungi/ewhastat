################################################################
#         R-codes for Sliced Inverse Regression
################################################################
################################################################
#          Power of matrix
################################################################
matpower=function(a,alpha){
  a = (a + t(a))/2
  tmp = eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%
           t(tmp$vectors))}

##############################################################
#      Center X (n by p matrix)        
##############################################################
center = function(x){
  return(t(t(x)-apply(x,2,mean)))}


################################################################
#           discretize
################################################################
discretize=function(y,h){
  n=length(y);m=round(n/h)
  y=y+.00001*mean(y)*rnorm(n)
  yord = y[order(y)]
  divpt=numeric();for(i in 1:(h-1)) divpt = c(divpt,yord[i*m+1])
  y1=rep(0,n);y1[y<divpt[1]]=1;y1[y>=divpt[h-1]]=h
  for(i in 2:(h-1)) y1[(y>=divpt[i-1])&(y<divpt[i])]=i
  return(y1)
}

################################################################
#          distance between subspaces
################################################################
dist=function(v1,v2){
  v1=as.matrix(v1);v2=as.matrix(v2)
  if(dim(v1)[1]>1){
    p1 <- v1%*%matpower(t(v1)%*%v1,-1)%*%t(v1)
    p2 <- v2%*%matpower(t(v2)%*%v2,-1)%*%t(v2)
    if(dim(v1)[1]==1){
      p1=v1%*%t(v1)/c(t(v1)%*%v1)
      p2=v2%*%t(v2)/c(t(v2)%*%v2)}
    d <- sqrt(sum((p1-p2)*(p1-p2)))}
  return(d)
}
################################################################
#                           sir
################################################################
sir=function(x,y,h,r,ytype){
  p=ncol(x);n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  ylabel=unique(ydis)
  prob=numeric();exy=numeric()
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n) 
  for(i in 1:h) exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))
  sirmat=t(exy)%*%diag(prob)%*%exy
  return(signrt%*%eigen(sirmat)$vectors[,1:r])}

################################################################
#                          Model 1
################################################################
n = 400
p = 10
sig = 0.02
x = matrix(rnorm(n*p),n,p)
y = x[,1]/(0.5+(x[,2]+1.5)^2) + sig*rnorm(n)
beta.true=matrix(0,p,2);beta.true[1,1]=beta.true[2,2]=1
##############################################################
#                         apply sir
##############################################################
h=10;r=2;ytype="continuous";beta=sir(x,y,h,r,ytype)
pred=center(x)%*%beta
dist(beta,beta.true)
beta



################################################################
#                          Model 2
################################################################
n = 100
p = 10
sig = 0.2
x = matrix(rnorm(n*p),n,p)
y = x[,1]*(x[,1]+x[,2]+1) + sig * rnorm(n)
beta.true=matrix(0,p,2);beta.true[1,1]=beta.true[2,2]=1
##############################################################
#                        apply sir
##############################################################
h=8;r=2;ytype="continuous";beta=sir(x,y,h,r,ytype)
pred=center(x)%*%beta
dist(beta,beta.true)
beta


################################################################
#                          Model 3
################################################################
n = 100
p = 10
sig = 0.2
x = matrix(rnorm(n*p),n,p)
y = x[,1] + x[,2] + x[,3] + x[,4] + sig * rnorm(n)
beta.true=matrix(0,p,1);beta.true[c(1,2,3,4),1]=1
##############################################################
#                         apply sir
##############################################################
h=8;r=1;ytype="continuous";beta=sir(x,y,h,r,ytype)
pred=center(x)%*%beta
beta

################################################################
#                          Model 4 (application)
################################################################
y = bigmac[,1]
x = bigmac[,2:10]
pairs(x)



############################################################## 
#              run sir on big mac data 
##############################################################
h=8;r=2;ytype="continuous";beta=sir(x,y,h,r,ytype)
par(mfrow=c(1,2))
pred = center(x) %*% beta 
plot(pred[,1],y)
plot(pred[,2],y)
par(mfrow=c(1,1))



############################################################## 
#              Penndigit data
##############################################################

xy.tra = pendigits.tra
xy.tra.069 = rbind(xy.tra[xy.tra[,17]==0,],
                   xy.tra[xy.tra[,17]==6,],xy.tra[xy.tra[,17]==9,])
x.tra=xy.tra.069[,1:16];y.tra=as.matrix(xy.tra.069[,17])
n=length(y.tra)

xy.tes = pendigits.tes
xy.tes.069 = rbind(xy.tes[xy.tes[,17]==0,],
                   xy.tes[xy.tes[,17]==6,],xy.tes[xy.tes[,17]==9,])
x.tes=xy.tes.069[,1:16];y.tes=xy.tes.069[,17]
n=length(y.tes)

##############################################################
#            training dataset
##############################################################
ytype="categorical";h=3; r=3
outcome <- sir(x.tra,y.tra,h,r,ytype)
pred <-  as.matrix(x.tra)%*%outcome
plot(pred[,1],pred[,2],xlab="first kernel PC",ylab="second kernel PC")
points(pred[y.tra==0,1],pred[y.tra==0,2],col="red")
points(pred[y.tra==6,1],pred[y.tra==6,2],col="green")
points(pred[y.tra==9,1],pred[y.tra==9,2],col="blue")


##############################################################
#            test dataset
##############################################################
pred <-  as.matrix(x.tes)%*%outcome
plot(pred[,1],pred[,2],xlab="first kernel PC",ylab="second kernel PC")
points(pred[y.tes==0,1],pred[y.tes==0,2],col="red")
points(pred[y.tes==6,1],pred[y.tes==6,2],col="green")
points(pred[y.tes==9,1],pred[y.tes==9,2],col="blue")



##############################################################
#            Another classification problem
##############################################################

library(e1071)
library(kernlab)
library(plot3D)

#
# 10-fold CV linear svm (using Kernel data)
#
svm.cvdim <- function(data, kdata, dim, k=10){
  m <- dim(data)[1]
  q <- dim(data)[2]
  if(dim==1){
    xdata <- as.matrix(kdata[,1])
  }
  else{
    xdata <- as.matrix(kdata[,1:dim])
  }
  y <- data[,1]  #was data[,1]
  cv.error.single <- c()
  n <- nrow(xdata)
  
  K <- k
  no.class <- length(unique(y))
  number <- table(y)
  #cat("number: ", number, "\n")
  index <- c()
  if(dim==1){
    xdata <- as.matrix(xdata[order(y),1])
  }
  else{
    xdata <- xdata[order(y),]
  }
  
  y <- sort(y)
  
  for(i in 1:no.class){
    samps <- sample(rep(1:K, length=number[i]), number[i], replace=FALSE)
    index <- c(index, samps)
  }
  ## stratified  CV
  for (i in 1:K){
    cv.x.test <- xdata[index==i,]
    cv.x.train <- xdata[index!=i,]
    cv.y.test <- y[index==i]
    cv.y.train <- y[index!=i]
    #svm.model <- svm(x=cv.x.train, y=as.factor(cv.y.train), cost=100, gamma=0.05)
    svm.model <- svm(x=cv.x.train, y=as.factor(cv.y.train), kernel="linear")
    predicted <- predict(svm.model, cv.x.test)
    svm.error <- mean(predicted != cv.y.test)
    cv.error.single <- c(cv.error.single, svm.error)
  }
  cv.error.m <- mean(cv.error.single)
  cat("cv.error.m: ", cv.error.m, "\n")
  cv.error.sd <- sd(cv.error.single)
  cat("cv.error.sd: ", cv.error.sd, "\n")
  return(list(cv.error.m,cv.error.sd))
}


#
# 10-fold CV linear svm (using original data)
#
svm.cv <- function(data, k=10, scale=FALSE){
  m <- dim(data)[1]
  q <- dim(data)[2]
  
  if(scale){
    xdata <- scale(data[,2:q])
  }
  else{
    xdata <- data[,2:q]
  }
  
  y <- data[,1]
  cv.error.single <- c()
  n <- nrow(xdata)
  
  K <- k
  no.class <- length(unique(y))
  number <- table(y)
  #cat("number: ", number, "\n")
  index <- c()
  xdata <- xdata[order(y),]
  y <- sort(y)
  
  for(i in 1:no.class){
    samps <- sample(rep(1:K, length=number[i]), number[i], replace=FALSE)
    index <- c(index, samps)
  }
  ## stratified  CV
  for (i in 1:K){
    cv.x.test <- xdata[index==i,]
    cv.x.train <- xdata[index!=i,]
    cv.y.test <- y[index==i]
    cv.y.train <- y[index!=i]
    #svm.model <- svm(x=cv.x.train, y=as.factor(cv.y.train), cost=100, gamma=0.05)
    svm.model <- svm(x=cv.x.train, y=as.factor(cv.y.train), kernel="linear")
    
    predicted <- predict(svm.model, cv.x.test)
    svm.error <- mean(predicted != cv.y.test)
    cv.error.single <- c(cv.error.single, svm.error)
  }
  cv.error.m <- mean(cv.error.single)
  cat("cv.error.m: ", cv.error.m, "\n")
  cv.error.sd <- sd(cv.error.single)
  cat("cv.error.sd: ", cv.error.sd, "\n")
  return(list(cv.error.m,cv.error.sd))
}

# gevd compute the generalized eigenvalue
# decomposition for ( a , b)
gevd <- function(a , b=diag(nrow(a))){
  bs <- mfunc(b , function(x) ginvx(sqrt(x)))
  ev <- eigen ( bs%*%a%*%bs )
  return(list(values=ev$values, vectors=bs%*%ev$vectors))
}
# ginvx i s a helper to compute reciprocals
ginvx <- function ( x ) { ifelse ( x==0,0,1/x ) }
# mfunc i s a helper to compute matrix functions
mfunc <- function ( a , fn=sqrt ) {
  e <- eigen(a);
  y <- e$vectors ; v <- e$values
  return(tcrossprod (y%*%diag(fn(v)),y))
}

parkinsons <- read.csv("~/Desktop/Desktop/Teaching/EWHA/다변량해석특론_2021_Fall/R codes/application/datasets/parkinsons.data", header=FALSE)
xdata <- scale(matrix(as.numeric(unlist(parkinsons[-1,-c(1,18)])),195,22))
ydata <- as.numeric(parkinsons[-1,18])
input.data <- cbind(ydata, xdata)



ytype="categorical";h=2; r=3
result <- sir(xdata,ydata,h,r,ytype)
pred <-  as.matrix(xdata)%*%result


plot(pred[,1], pred[,2], xlab='first predictor', ylab='second predictor')
points(pred[ydata==0,1],pred[ydata==0,2],col="red",pch=1)
points(pred[ydata==1,1],pred[ydata==1,2],col="blue",pch=1)

no.dim <- 2
svm.cvdim(input.data, pred, dim=no.dim)
svm.cv(input.data,k=10,FALSE) 





