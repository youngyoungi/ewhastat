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
  }
    if(dim(v1)[1]==1){
      p1=v1%*%t(v1)/c(t(v1)%*%v1)
      p2=v2%*%t(v2)/c(t(v2)%*%v2)}
    d <- sqrt(sum((p1-p2)*(p1-p2)))
  return(d)
}

##############################################################
#                   symmtrize a matrix
############################################################## 
symmetry = function(a){
  return((a + t(a))/2)}


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
#                          pir 
################################################################
pir=function(x,y,m,r){
  xc=t(t(x)-apply(x,2,mean))
  signrt=matpower(var(x),-1/2)
  xstand = xc%*%signrt
  f=numeric();ystand=(y-mean(y))/sd(y)
  for(i in 1:m) f = cbind(f, ystand^i) 
  sigxf=cov(xstand,f);sigff=var(f) 
  cand=sigxf%*%solve(sigff)%*%t(sigxf)
  return(signrt%*%eigen(symmetry(cand))$vectors[,1:r])
}
################################################################
#                          save
################################################################
save=function(x,y,h,r,ytype){
  p=ncol(x);n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  ylabel=unique(ydis)
  prob=numeric() 
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy = array(0,c(p,p,h))
  for(i in 1:h) vxy[,,i] = var(xst[ydis==ylabel[i],]) 
  savemat=0
  for(i in 1:h){
    savemat=savemat+prob[i]*(vxy[,,i]-diag(p))%*%(vxy[,,i]-diag(p))}
  return(signrt%*%eigen(savemat)$vectors[,1:r])}


################################################################
#                        sir ii 
################################################################
sirii=function(x,y,h,r,ytype="continuous"){                                         
  p=ncol(x);n=nrow(x)                                              
  signrt=matpower(var(x),-1/2)                                     
  xst=t(t(x)-apply(x,2,mean))%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)  
  if(ytype=="categorical") ydis=y 
  ylabel=unique(ydis)
  prob=numeric();exy=numeric()                                     
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)               
  for(i in 1:h) exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))         
  sirmat=t(exy)%*%diag(prob)%*%exy                                 
  vxy = array(0,c(p,p,h))                                          
  for(i in 1:h) vxy[,,i] = var(xst[ydis==ylabel[i],])                      
  savemat=0                                                        
  for(i in 1:h){                                                   
    savemat=savemat+prob[i]*(vxy[,,i]-diag(p))%*%(vxy[,,i]-diag(p))} 
  siriimat=savemat-sirmat%*%t(sirmat)                              
  return(signrt%*%eigen(siriimat)$vectors[,1:r])}
################################################################
#                               cr 
################################################################
cr=function(x,y,percent,r){                                  
  tradeindex12 = function(k,n){                                
    j = ceiling(k/n)                                             
    i = k - (j-1)*n                                              
    return(c(i,j))}                                                                                    
  mu=apply(x,2,mean);signrt=matpower(var(x),-1/2)              
  z=t(t(x)-mu)%*%signrt                                        
  n=dim(x)[1];p = dim(x)[2]                                    
  ymat=matrix(y,n,n)                                           
  deltay=c(abs(ymat - t(ymat)))                                
  singleindex=(1:n^2)[deltay < percent*mean(deltay)]           
  contourmat=matrix(0,p,p)                                     
  for(k in singleindex){                                       
    doubleindex=tradeindex12(k,n)                                
    deltaz=z[doubleindex[1],]-z[doubleindex[2],]                 
    contourmat=contourmat+deltaz %*% t(deltaz)}                  
  signrt=matpower(var(x),-1/2)                                 
  return(signrt%*%eigen(contourmat)$vectors[,p:(p-r+1)])                   
}
################################################################
#                                 dr 
################################################################
dr=function(x,y,h,r,ytype){
  p=ncol(x);n=nrow(x)
  signrt=matpower(var(x),-1/2)
  xc=t(t(x)-apply(x,2,mean))
  xst=xc%*%signrt
  if(ytype=="continuous") ydis=discretize(y,h)
  if(ytype=="categorical") ydis=y
  ylabel=unique(ydis)
  prob=numeric() 
  for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy = array(0,c(p,p,h));exy=numeric()
  for(i in 1:h) {
    vxy[,,i]=var(xst[ydis==ylabel[i],])
    exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))}
  mat1 = matrix(0,p,p);mat2 = matrix(0,p,p)
  for(i in 1:h){
    mat1 = mat1+prob[i]*(vxy[,,i]+exy[i,]%*%t(exy[i,]))%*%
      (vxy[,,i]+exy[i,]%*%t(exy[i,]))
    mat2 = mat2+prob[i]*exy[i,]%*%t(exy[i,])}
  out = 2*mat1+2*mat2%*%mat2+2*sum(diag(mat2))*mat2-2*diag(p)
  return(signrt%*%eigen(out)$vectors[,1:r])
}

#####################################################
#                       OLS              
##################################################### 
ols = function(x,y){
  return(solve(var(x))%*%cov(x,y))}
beta=ols(x,y)
pred=c(t(t(x)-apply(x,2,mean))%*%beta)
plot(pred,y,xlab="OLS predictor",ylab="y")
cor(rank(y),rank(pred))
beta1=ols(xlam,y)
pred1=c(t(t(xlam)-apply(xlam,2,mean))%*%beta1)
plot(pred1,y,xlab="OLS predictor",ylab="y")
cor(rank(y),rank(pred1))


#####################################################
#                     PhD             
##################################################### 
phd = function(x,y,r){
  matpower=function(a,alpha){
    a=(a+t(a))/2;tmp=eigen(a)
    return(tmp$vectors%*%diag((tmp$values)^alpha)%*%t(tmp$vectors))}
  n=length(y);signrt=matpower(var(x),-1/2)
  z = center(x)%*%signrt
  return(signrt%*%eigen(t(z*(y-mean(y)))%*%z/n)$vectors[,1:r])}

#####################################################
#                    e-PhD            
##################################################### 
phd1 = function(x,y,r){
  n = length(y)
  sig = var(x)
  signrt = matpower(sig,-1/2)
  z = center(x)%*%signrt
  yc = y - mean(y) - center(x)%*%solve(sig)%*%cov(x,y)
  return(signrt%*%eigen(t(z*y)%*%z/n)$vectors[,1:r])
}


##########################################################
#                     iht
##########################################################
iht=function(x,y,r){
  standmat=function(x){
    mu=apply(x,2,mean);sig=var(x);signrt=matpower(sig,-1/2)
    return(t(t(x)-mu)%*%signrt)}
  z=standmat(x);szy=cov(z,y)
  p=dim(z)[2];n=dim(z)[1];szzy=t(z)%*%(z*y)/n;imat=szy
  for(i in 1:(p-1)) imat=cbind(imat,matpower(szzy,i)%*%szy)
  sxx=var(x);sqrtsxx=matpower(sxx,-1/2)
  return(sqrtsxx%*%eigen(imat%*%t(imat))$vectors[,1:r])}
  
  



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


############################################################## 
#                        apply pir 
##############################################################
m=6;r=2;ytype="continuous";beta=pir(x,y,m,r)
beta

############################################################## 
#                        run ols
##############################################################
m=1;r=2;ytype="continuous";beta=pir(x,y,m,r)
beta

##############################################################
#                        apply save
##############################################################
h=2;r=2;ytype="continuous";beta=save(x,y,h,r,ytype)
beta


##############################################################
#                         apply sirii
##############################################################
h=10;r=2;ytype="continuous";beta=sirii(x,y,h,r,ytype)
pred=center(x)%*%beta
dist(beta,beta.true)
beta


##############################################################
#                         apply cr
##############################################################
beta=cr(x,y,0.01,r)
pred=center(x)%*%beta
dist(beta,beta.true)
beta

##############################################################
#                         apply dr
##############################################################
h=10;r=2;ytype="continuous";beta=dr(x,y,h,r,ytype)
pred=center(x)%*%beta
dist(beta,beta.true)
beta

##############################################################
#                         apply ols
##############################################################
beta=ols(x,y)
pred=center(x)%*%beta
dist(beta,beta.true)
beta

##############################################################
#                         apply phd
##############################################################
beta=phd(x,y,2)
pred=center(x)%*%beta
dist(beta,beta.true)
beta


##############################################################
#                         apply iht
##############################################################
beta=iht(x,y,2);pred=center(x)%*%beta
plot(pred[,1],y,xlab="first IHT direction",ylab="y")
plot(pred[,2],y,xlab="second IHT direction",ylab="y")

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

############################################################## 
#                        apply pir 
##############################################################
m=6;r=2;ytype="continuous";beta=pir(x,y,m,r)
beta

############################################################## 
#                        run ols
##############################################################
m=1;r=2;ytype="continuous";beta=pir(x,y,m,r)
beta

##############################################################
#                        apply save
##############################################################
h=2;r=2;ytype="continuous";beta=save(x,y,h,r,ytype)
beta

##############################################################
#                         apply sirii
##############################################################
h=10;r=2;ytype="continuous";beta=sirii(x,y,h,r,ytype)
pred=center(x)%*%beta
dist(beta,beta.true)
beta


##############################################################
#                         apply cr
##############################################################
beta=cr(x,y,0.01,r)
pred=center(x)%*%beta
dist(beta,beta.true)
beta

##############################################################
#                         apply dr
##############################################################
h=10;r=2;ytype="continuous";beta=dr(x,y,h,r,ytype)
pred=center(x)%*%beta
dist(beta,beta.true)
beta


##############################################################
#                         apply ols
##############################################################
beta=ols(x,y)
pred=center(x)%*%beta
dist(beta,beta.true)
beta


##############################################################
#                         apply phd
##############################################################
beta=phd(x,y,2)
pred=center(x)%*%beta
dist(beta,beta.true)
beta

##############################################################
#                         apply iht
##############################################################
beta=iht(x,y,2);pred=center(x)%*%beta
plot(pred[,1],y,xlab="first IHT direction",ylab="y")
plot(pred[,2],y,xlab="second IHT direction",ylab="y")




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

############################################################## 
#                        apply pir 
##############################################################
m=6;r=2;ytype="continuous";beta=pir(x,y,m,r)
beta

############################################################## 
#                        run ols
##############################################################
m=1;r=2;ytype="continuous";beta=pir(x,y,m,r)
beta

##############################################################
#                        apply save
##############################################################
h=2;r=2;ytype="continuous";beta=save(x,y,h,r,ytype)
beta


##############################################################
#                         apply sirii
##############################################################
h=10;r=1;ytype="continuous";beta=sirii(x,y,h,r,ytype)
pred=center(x)%*%beta
beta


##############################################################
#                         apply cr
##############################################################
beta=cr(x,y,0.01,r)
pred=center(x)%*%beta
beta

##############################################################
#                         apply dr
##############################################################
h=10;r=2;ytype="continuous";beta=dr(x,y,h,r,ytype)
pred=center(x)%*%beta
beta

##############################################################
#                         apply ols
##############################################################
beta=ols(x,y)
pred=center(x)%*%beta
dist(beta,beta.true)
beta


##############################################################
#                         apply phd
##############################################################
beta=phd(x,y,2)
pred=center(x)%*%beta
dist(beta,beta.true)
beta

##############################################################
#                         apply iht
##############################################################
beta=iht(x,y,2);pred=center(x)%*%beta
plot(pred[,1],y,xlab="first IHT direction",ylab="y")
plot(pred[,2],y,xlab="second IHT direction",ylab="y")



################################################################
#                          Model 4 (symmetric)
################################################################
n = 100
p = 10
sig = 0.2
x = matrix(rnorm(n*p),n,p)
y = x[,1]^2 + sig * rnorm(n)
beta.true=matrix(0,p,1);beta.true[c(1,2,3,4),1]=1
##############################################################
#                         apply sir
##############################################################
h=8;r=1;ytype="continuous";beta=sir(x,y,h,r,ytype)
pred=center(x)%*%beta
beta

############################################################## 
#                        apply pir 
##############################################################
m=6;r=1;ytype="continuous";beta=pir(x,y,m,r)
beta

############################################################## 
#                        run ols
##############################################################
m=1;r=1;ytype="continuous";beta=pir(x,y,m,r)
beta

##############################################################
#                        apply save
##############################################################
h=2;r=1;ytype="continuous";beta=save(x,y,h,r,ytype)
beta

##############################################################
#                         apply sirii
##############################################################
h=10;r=2;ytype="continuous";beta=sirii(x,y,h,r,ytype)
pred=center(x)%*%beta
beta


##############################################################
#                         apply cr
##############################################################
beta=cr(x,y,0.01,r)
pred=center(x)%*%beta
beta

##############################################################
#                         apply dr
##############################################################
h=10;r=2;ytype="continuous";beta=dr(x,y,h,r,ytype)
pred=center(x)%*%beta
beta

##############################################################
#                         apply ols
##############################################################
beta=ols(x,y)
pred=center(x)%*%beta
dist(beta,beta.true)
beta


##############################################################
#                         apply phd
##############################################################
beta=phd(x,y,2)
pred=center(x)%*%beta
dist(beta,beta.true)
beta

##############################################################
#                         apply iht
##############################################################
beta=iht(x,y,2);pred=center(x)%*%beta
plot(pred[,1],y,xlab="first IHT direction",ylab="y")
plot(pred[,2],y,xlab="second IHT direction",ylab="y")


############################################################## 
#         Compare with SIR, SAVE, SIRII, and CR on 
#        for situations where both a asymmetric trend and 
#        and a symmetric trend are present 
##############################################################
##############################################################
#         generate model 
##############################################################
n = 100
p = 10
sig = 0.2
x = matrix(rnorm(n*p),n,p)
y = x[,1]^2 + 2 * sin(x[,2]) + sig * rnorm(n)
beta.true=matrix(0,p,2);beta.true[1,1]=beta.true[2,2]=1
spin(y,x[,1:2],900,2)
##############################################################
#      apply sir
##############################################################
h=8;r=2;ytype="continuous";beta=sir(x,y,h,r,ytype)
pred=center(x)%*%beta
spin(y,pred,500,3)
dist(beta,beta.true)
##############################################################
#       apply save
##############################################################
h=2;r=2;ytype="continuous";beta=save(x,y,h,r,ytype)
pred=center(x)%*%beta
spin(y,pred,900,3)
pairs(cbind(pred,y))
dist(beta,beta.true)
##############################################################
#        apply sirii
##############################################################
h=2;r=2;ytype="continuous";beta=sirii(x,y,h,r,ytype)
pred=center(x)%*%beta
spin(y,pred,200,3)
pairs(cbind(pred,y))
dist(beta,beta.true)
##############################################################
#        apply cr
##############################################################
percent=0.15;r=2;beta=cr(x,y,percent,r)
pred=center(x)%*%beta
spin(y,pred,200,3)
pairs(cbind(pred,y))
dist(beta,beta.true)
##############################################################
#         apply dr
##############################################################
h=4;r=2;ytype="continuous";beta=dr(x,y,h,r,ytype)
pred=center(x)%*%beta
spin(y,pred,200,3)
pairs(cbind(pred,y))
dist(beta,beta.true)


##############################################################
#                         apply ols
##############################################################
beta=ols(x,y)
pred=center(x)%*%beta
dist(beta,beta.true)
beta

##############################################################
#                         apply phd
##############################################################
beta=phd(x,y,2)
pred=center(x)%*%beta
dist(beta,beta.true)
beta

##############################################################
#                         apply iht
##############################################################
beta=iht(x,y,2);pred=center(x)%*%beta
plot(pred[,1],y,xlab="first IHT direction",ylab="y")
plot(pred[,2],y,xlab="second IHT direction",ylab="y")






############################################################## 
#  demonstration of empirical directions 
##############################################################
############################################################## 
#  generate data
##############################################################
n=50;p=2;x=rnorm(n*p);x=matrix(x,n,p);eps=rnorm(n)
y=sin(x[,1])+0.3*eps
plot(x[,1],y)
############################################################## 
# 11.2 empirical directions
##############################################################
plot(x[,1],x[,2])
for(i in 2:n) for(j in 1:(i-1)) lines(c(x[i,1],x[j,1]),c(x[i,2],x[j,2]))
############################################################## 
# 11.3 contour directions
##############################################################
plot(x[,1],x[,2])
for(i in 2:n) for(j in 1:(i-1)) {
  if(abs(y[i]-y[j])<0.1) lines(c(x[i,1],x[j,1]),c(x[i,2],x[j,2]))}





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
#              run pir on big mac data 
##############################################################
m=4;r=2;ytype="continuous";beta=pir(x,y,m,r)
par(mfrow=c(1,2))
pred = center(x) %*% beta 
plot(pred[,1],y)
plot(pred[,2],y)
par(mfrow=c(1,1))


############################################################## 
#              run save on big mac data 
##############################################################
h=8;r=2;ytype="continuous";beta=save(x,y,h,r,ytype)
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
#            training dataset (sir)
##############################################################
ytype="categorical";h=3; r=3
outcome <- sir(x.tra,y.tra,h,r,ytype)
pred <-  as.matrix(x.tra)%*%outcome
plot(pred[,1],pred[,2],xlab="first kernel PC",ylab="second kernel PC")
points(pred[y.tra==0,1],pred[y.tra==0,2],col="red")
points(pred[y.tra==6,1],pred[y.tra==6,2],col="green")
points(pred[y.tra==9,1],pred[y.tra==9,2],col="blue")


##############################################################
#            test dataset  (sir)
##############################################################
pred <-  as.matrix(x.tes)%*%outcome
plot(pred[,1],pred[,2],xlab="first kernel PC",ylab="second kernel PC")
points(pred[y.tes==0,1],pred[y.tes==0,2],col="red")
points(pred[y.tes==6,1],pred[y.tes==6,2],col="green")
points(pred[y.tes==9,1],pred[y.tes==9,2],col="blue")


##############################################################
#            training dataset (pir)
##############################################################
m=2; r=3
outcome <- pir(x.tra,y.tra,m,r)
pred <-  as.matrix(x.tra)%*%outcome
plot(pred[,1],pred[,2],xlab="first kernel PC",ylab="second kernel PC")
points(pred[y.tra==0,1],pred[y.tra==0,2],col="red")
points(pred[y.tra==6,1],pred[y.tra==6,2],col="green")
points(pred[y.tra==9,1],pred[y.tra==9,2],col="blue")


##############################################################
#            test dataset  (pir)
##############################################################
pred <-  as.matrix(x.tes)%*%outcome
plot(pred[,1],pred[,2],xlab="first kernel PC",ylab="second kernel PC")
points(pred[y.tes==0,1],pred[y.tes==0,2],col="red")
points(pred[y.tes==6,1],pred[y.tes==6,2],col="green")
points(pred[y.tes==9,1],pred[y.tes==9,2],col="blue")

##############################################################
#            training dataset (save)
##############################################################
ytype="categorical";h=3; r=3
outcome <- save(x.tra,y.tra,h,r,ytype)
pred <-  as.matrix(x.tra)%*%outcome
plot(pred[,1],pred[,2],xlab="first kernel PC",ylab="second kernel PC")
points(pred[y.tra==0,1],pred[y.tra==0,2],col="red")
points(pred[y.tra==6,1],pred[y.tra==6,2],col="green")
points(pred[y.tra==9,1],pred[y.tra==9,2],col="blue")


##############################################################
#            test dataset  (save)
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

parkinsons <- read.csv("~/Desktop/sdr class/sir/parkinsons.data", header=FALSE)
xdata <- scale(matrix(as.numeric(unlist(parkinsons[-1,-c(1,18)])),195,22))
ydata <- as.numeric(parkinsons[-1,18])
input.data <- cbind(ydata, xdata)

##############################################################
#            training dataset (sir)
##############################################################

ytype="categorical";h=2; r=3
result <- sir(xdata,ydata,h,r,ytype)
pred <-  as.matrix(xdata)%*%result


plot(pred[,1], pred[,2], xlab='first predictor', ylab='second predictor')
points(pred[ydata==0,1],pred[ydata==0,2],col="red",pch=1)
points(pred[ydata==1,1],pred[ydata==1,2],col="blue",pch=1)



no.dim <- 2
svm.cvdim(input.data, pred, dim=no.dim)
svm.cv(input.data,k=10,FALSE)  #naive case

##############################################################
#            training dataset (save)
##############################################################

ytype="categorical";h=2; r=3
result <- save(xdata,ydata,h,r,ytype)
pred <-  as.matrix(xdata)%*%result


plot(pred[,1], pred[,2], xlab='first predictor', ylab='second predictor')
points(pred[ydata==0,1],pred[ydata==0,2],col="red",pch=1)
points(pred[ydata==1,1],pred[ydata==1,2],col="blue",pch=1)

no.dim <- 2
svm.cvdim(input.data, pred, dim=no.dim)
svm.cv(input.data,k=10,FALSE) 





