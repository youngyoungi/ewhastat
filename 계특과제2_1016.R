# 준비
library(numDeriv)


f_simple<-function(x){
  cos(x)
}

f1<-function(x){
  log(x)/(1+x)
}

f_dif<-function(x){
  (x^3+1)/log(1/x)
}

# NR-method 구현

NR<-function(f, x0, thres=10^(-10),N=1000){
  k<-0  #iteration
  err<-1
     while (k<=N && err>=thres){   #stopping rule: err가 thres보다 작아지면 멈춘다.
      z<-genD(f,x0)
      h_x<- (-z$D[1,1]/z$D[1,2])
      x<-x0+h_x
      oldx0<-x0  #for calculating err
      x0<-x  #initial point update
      err<-abs(x0-oldx0) #err update
      k<-k+1 #iteration
  }
  res=genD(f,x0)
  print(paste("number of iteration=",k,sep=''))
  print(paste("result(x0)=",x0,sep=""))
  print(paste("g'(x0)=",res$D[1,1],sep=''))
  print(paste("g(x0)=",res$f0,sep=""))
  
}


system.time(NR(f_simple,4))
system.time(NR(f1,2))
system.time(NR(f_dif,2))

# secant method 구현

secant<-function(f, x1,x0, thres=10^(-10),N=500){
  #x0<-max(x0,x1)
  #x1<-min(x0,x1)
  z<-genD(f,c(x0,x1))
  h_x<- (-z$D[2,2]*(x1-x0)/(z$D[2,2]-z$D[1,1]))
  x<-x1+h_x
  k<-0  #iteration
  err<-1
  while (k<=N && err>=thres){   #stopping rule: err가 thres보다 작아지면 멈춘다.
    z0<-genD(f,x0)
    z1<-genD(f,x1)
    
    h_x<- (-z1$D[1]*(x1-x0)/(z1$D[1]-z0$D[1]))
    x<-x1+h_x
    x0<-x1
    x1<-x
    oldx1<-x0
    err<-abs(x0-x1) #err update
    k<-k+1 #iteration
  }
  res=genD(f,x1)
  print(paste("number of iteration=",k,sep=''))
  print(paste("result(x*)=",x1,sep=""))
  print(paste("g'(x)=",res$D[1,1],sep=''))
  print(paste("g(x)=",res$f0,sep=""))
  
}

system.time(secant(f_simple,1,4))
system.time(secant(f1,2,4))
system.time(secant(f_dif,1,10))
