source("main.R")

#### parameter setting
gamma=4 # kernel smoothing bandwidth
lam1 = 1/(2*gamma^2)
lam2 = 3/(4*gamma^4)

##############################################################
####  Checking excursion length distribution
#############################################################
## f_u(l): Eq 5.6
fx = function(x,lam1,lam2,u){
  if(lam1<=0 | lam2<=0) stop("lambda1 and lambda2 must be postive numbers")
  det = lam2-lam1^2
  den = lam2*x^4-16*lam1*x^2+64
  if(det<=0) stop("lambda2 is smaller than lambda1 square")
  num = lam2*x^2-8*lam1
  p1 = sqrt(lam2/(2*pi))*(1-pnorm(sqrt(lam2/det)*u)) +
    lam1*dnorm(u)*(1-pnorm(-lam1*u/sqrt(det)))
  p2 = 128*x*det/den^(3/2)*dnorm(8*u/sqrt(den)) # correct version
  p3 = (1+u^2*num^2/(den*det))*(1-pnorm(num*u/(sqrt(det*den))))
  p4 = u*num/(sqrt(det*den))*dnorm(u*num/sqrt(det*den))
  return(1/p1*p2*(p3-p4))
}
## F_u(l): Eq 5.5
Fx = function(x,lam1,lam2,u){
  if(lam1<=0 | lam2<=0) stop("lambda1 and lambda2 must be postive numbers")
  det = lam2-lam1^2
  den = lam2*x^4-16*lam1*x^2+64
  if(det<=0) stop("lambda2 is smaller than lambda1 square")
  num = lam2*x^2-8*lam1
  p1 = sqrt(lam2/(2*pi))*(1-pnorm(sqrt(lam2/det)*u))
  p2 = num/sqrt(den)*dnorm(8*u/sqrt(den))*(1-pnorm(num*u/sqrt(det*den)))
  p3 = lam1*dnorm(u)*(1-pnorm(-lam1*u/sqrt(det)))
  return((p1-p2)/(p1+p3))
}
## verification 
size=5000000; gamma=4
data = rnorm(size)
sdata = smth.gau(data,gamma); sdata = sdata/sd(sdata) # standardize
locmax = which.peaks(sdata)
par(mfrow=c(2,2),mar=c(3,2,2,1))
for (u in c(2.5,3,3.2,3.5)){
  ulocmax = locmax[sdata[locmax]>u]
  width = sapply(lapply(ulocmax,search,u=u,data=sdata),sum)
  ## tail probability
  plot(x =ecdf.tail(width)[,1],y = ecdf.tail(width)[,2],type="l",lwd=2,
       xlab="Excursion Length",ylab="Tail Probability",main=paste("u = ",u,sep=""))
  x = seq(min(width),max(width),length.out=100)
  lines(x = x,y=Fx(x,lam1,lam2,u),col="red",lwd=2)
  legend("topright",legend=c("Numerical tail","Theoretical tail"),
         col=c("black","red"),lwd=c(2,2),bty="n")
}
#dev.off()
#############################################################
#### checking excursion area distribution
#############################################################
gv = function(v,lam1,lam2,u){
  if(lam1<=0 | lam2<=0) stop("lambda1 and lambda2 must be postive numbers")
  det = lam2-lam1^2
  p1 = sqrt(lam2/(2*pi))*(1-pnorm(sqrt(lam2/det)*u))
  p2 = lam1*dnorm(u)*(1-pnorm(-lam1*u/sqrt(det)))
  p3 = 2048/(81*sqrt(2*pi*det)*v^5)
  f = function(x) x^6*exp(-(lam1*(x+u)-32*x^3/(9*v^2))^2/(2*det))*dnorm(x+u)
  p4 = integrate(f=f,lower=0,upper=Inf)$value
  # integrate vs. trapz
  return(1/(p1+p2)*p3*p4)
}
## 
Gv = function(V,lam1,lam2,u){
  X = seq(u,20,0.1)
  n = 10000
  Z = rnorm(n,0,1)
  f1 = function(x){
    mZ = abs(sqrt(lam2-lam1^2)*Z-lam1*x)
    ind1 = which(Z>=x*lam1/sqrt(lam2-lam1^2))
    mZ[ind1] = 0
    return(mean(mZ*dnorm(x)))
  }
  f2 = function(x,V){
    mZ = abs(sqrt(lam2-lam1^2)*Z-lam1*x)
    ind1 = which(Z>=x*lam1/sqrt(lam2-lam1^2))
    ind2 = which(Z<=(x*lam1 - 32*(x-u)^3/(9*V^2))/sqrt(lam2-lam1^2))
    mZ[c(ind1,ind2)] = 0
    return(mean(mZ)*dnorm(x))
  }
  temp1 = trapz(x=X,y=sapply(X,f1))
  temp2 = trapz(x=X,y=sapply(X,f2,V=V))
  return(temp2/temp1)
}
## verification
locmax = which.peaks(sdata)
par(mfrow=c(2,2),mar=c(3,2,2,1))
for (u in c(2.5,3,3.2,3.5)) {
  ulocmax = locmax[sdata[locmax]>u]
  area = sapply(ulocmax,auc,u=u,data=sdata); area = area[!is.na(area)]
  vc = seq(0,20,0.5)
  vGv = sapply(vc,Gv,lam1=lam1,lam2=lam2,u=u)
  plot(vc,vGv,xlab="Excursion Area",ylab="Tail Probability",type="l",lwd=2, col="red",
       main=paste("u = ",u,sep = ""))
  abline(h=0,lty=2,lwd=1,col="grey")
  lines(x = ecdf.tail(area)[,1],y=ecdf.tail(area)[,2],lwd=2)
  
  legend("topright",legend=c("Theoretical tail","Numerical tail"),lwd=c(2,2),
         col=c("red","black"),bty="n")
}

###############################################################
#### 2D: checking excursion extent distribution
###############################################################
F2x = function(l,u,n=1000,gamma=4){
  # denominator
  F1 = function(x){
    Lam = matrix(rnorm(n*2),ncol=2)
    ind  = (Lam[,1] < Lam[,2] & Lam[,2] < x/sqrt(2)) + 0
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2)) * ind) * dnorm(x)
  }
  # numerator
  F2 = function(x,l){
    Lam = matrix(rnorm(n*2),ncol=2)
    ind1 = (Lam[,1] < Lam[,2] & Lam[,2] < x/sqrt(2)) + 0
    ind2 = (4*pi^2*(x-u)^2 >= l^2*(Lam[,2]-x/sqrt(2))*(Lam[,1]-x/sqrt(2))/(2*gamma^4)) + 0
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2))*ind1*ind2)*dnorm(x)
  }
  vx = seq(u,10,0.1) 
  temp2 = trapz(vx,sapply(vx,F2,l=l))
  temp1 = trapz(vx,sapply(vx,F1))
  return(temp2/temp1)
}
##
nr = 5000; nc = 5000   # 10000 * 10000
noise = matrix(rnorm(nr*nc),nr,nc)
sdata = smoothie::kernel2dsmooth(noise,"gauss",sigma=4,nx=nr,ny=nc)
sdata = sdata/0.0704  # 0.0704
# sigma = sqrt(1/(4*pi*gamma^2)) 
for (u in c(3,3.2,3.5,3.8)) {
  xlim = 0:50
  theoval = sapply(xlim,F2x,u=u)
  plot(predict(loess(theoval~xlim)),lwd=2, type="l",ylab="Tail Probability", col="red",
       xlab="excursion area",main=paste(" u = ",u,sep=""))
  clusters = cross.area.cluster(sdata,u)
  lines(x = sort(clusters$area_lower$area),y=ecdf.tail(clusters$area_lower$area)[,2],lwd=2,col="blue")
  lines(x = sort(clusters$area_fit$area),y=ecdf.tail(clusters$area_fit$area)[,2],lwd=2,col="black")
  abline(h=0.1,lwd=1,col="grey")
  legend("topright",legend=c("theoretical","lower","fitted"),
         lwd=rep(2,3),col=c("red","blue","black"),bty="n")
}  

###############################################################
#### 2D: checking excursion volume distribution
###############################################################
G2v = function(l,u,n=1000,gamma=4){
  # denominator
  F1 = function(x){
    Lam = matrix(rnorm(n*2),ncol=2)
    ind  = (Lam[,1] < Lam[,2] & Lam[,2] < x/sqrt(2)) + 0
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2)) * ind) * dnorm(x)
  }
  # numerator
  F2 = function(x,l){
    Lam = matrix(rnorm(n*2),ncol=2)
    ind1 = (Lam[,1] < Lam[,2] & Lam[,2] < x/sqrt(2)) + 0
    ind2 = (pi^2*(x-u)^4 >= l^2*(Lam[,2]-x/sqrt(2))*(Lam[,1]-x/sqrt(2))/(2*gamma^4)) + 0
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2))*ind1*ind2)*dnorm(x)
  }
  vx = seq(u,10,0.1) 
  temp2 = trapz(vx,sapply(vx,F2,l=l))
  temp1 = trapz(vx,sapply(vx,F1))
  return(temp2/temp1)
}
##
for (u in c(3,3.2,3.5,3.8)) {
  xlim = 0:40
  theoval = sapply(xlim,G2v,u=u)
  plot(x=xlim,y=theoval,lwd=2, type="l",ylab="Tail Probability", col="red",
       xlab="excursion area",main=paste(" u = ",u,sep=""))
  xx = cross.area.cluster(sdata,u)
  inner_fill=do.call("rbind",lapply(split(xx$cluster,xx$cluster$group),fill.point,data = sdata))
  inner_vol = do.call("rbind",lapply(split(inner_fill,inner_fill$group),getVolume)); inner_vol = inner_vol[!is.na(inner_vol$volume),"volume"]
  fit_fill = do.call("rbind",lapply(split(xx$fit_point,xx$fit_point$group),fill.point,data = sdata,isfit=T,clusters=xx$cluster))
  fit_vol = do.call("rbind",lapply(split(fit_fill,fit_fill$group),getVolume)); fit_vol = fit_vol[!is.na(fit_vol$volume),"volume"]
  #
  lines(x = sort(inner_vol),y=ecdf.tail(inner_vol)[,2],lwd=2,col="blue")
  lines(x = sort(fit_vol),y=ecdf.tail(fit_vol)[,2],lwd=2)
  abline(h=0.1,lwd=1,col="grey")
  legend("topright",legend=c("theoretical","lower","fitted"),
         lwd=rep(2,3),col=c("red","blue","black"),bty="n")
}  

###############################################################
#### Joint Distribution
###############################################################