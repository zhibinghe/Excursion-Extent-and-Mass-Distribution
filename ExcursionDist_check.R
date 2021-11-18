###### Verify some theoretical results of Excursion Distribution
library(dSTEM)   ## install from github  zhibinghe/ChangePoint
library(pracma)  ## compute area under a curve
#### Basic functions
## compute tail probability based on empirical cdf
ecdf.tail = function(x){
  fecdf = ecdf(x)
  x = sort(x)
  cbind(x,1 - fecdf(x))
}
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
## search endpoints of f(x)=u
search = function(x,u,data){
  # x: locatio of the local maximum 
  # left search
  lcount = 1
  while(data[x-lcount]>u & lcount < x-1){ lcount = lcount + 1}
  lcount = ifelse(lcount==x-1,NA,lcount)
  lcount = lcount - 1 + (sdata[x-lcount+1]-u)/(data[x-lcount+1]-data[x-lcount])
  # right search
  rcount = 1
  while(data[x+rcount]>u & rcount < length(data)-x){ rcount = rcount + 1}  
  rcount = ifelse(rcount==length(data)-x,NA,rcount)
  rcount = rcount -1 + (data[x+rcount-1]-u)/(data[x+rcount-1]-data[x+rcount])
  return(c(lcount,rcount))
}
#### Theoretical distribution
u=3
gamma=4
lam1 = 1/(2*gamma^2)
lam2 = 3/(4*gamma^4)
x = seq(0.5,10,by=0.1)
plot(x=x,y= 1 - Fx(x,lam1,lam2,u),ylim = c(0,1),type="l",lwd=2,main="CDF")
abline(h=1,col="red",lty=2)
## f(x) vs. dF(x)
plot(x,fx(x,lam1,lam2,u),lwd=2,type="l",xlab="Excursion length")
lines(x,y=-numDeriv::grad(Fx,x,lam1=lam1,lam2=lam2,u=u),lty=2,col="red")
legend("topright",legend=c("dF'(x)","f(x)"),col=c("red","black"),lty=c(2,1),bty="n")
#### Numerical distribution
size=10000000; gamma=4
data = rnorm(size)
sdata = smth.gau(data,gamma); sdata = sdata/sd(sdata) # standardize
locmax = which.peaks(sdata)
ulocmax = locmax[sdata[locmax]>u]
width = sapply(lapply(ulocmax,search,u=u,data=sdata),sum)
## tail probability
plot(x =ecdf.tail(width)[,1],y = ecdf.tail(width)[,2],type="l",lwd=2,
     xlab="Excursion Length",ylab="Tail Probability",main=paste("u = ",u,sep=""))
x = seq(min(width),max(width),length.out=100)
lines(x = x,y=Fx(x,lam1,lam2,u),col="red",lwd=2)
legend("topright",legend=c("Numerical tail","Theoretical F(x)"),
       col=c("black","red"),lwd=c(2,2),bty="n")
#### The distribution of Excursion length is correct!!!
#############################################################
#### Checking excursion area distribution
#############################################################
#### Theoretical distribution
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
#### Numerical Distribution
## using the same data as the previous one
auc = function(x,u,data){
  # x: location of local maxima
  endpoints = search(x,u,data)
  temp = (x-round(endpoints[1])):(x+round(endpoints[2]))
  traparea = trapz(x=temp,y=sdata[temp])- u*(tail(temp,1)-temp[1])
  traparea = ifelse(is.na(traparea),0,traparea)
}
locmax = which.peaks(sdata)
ulocmax = locmax[sdata[locmax]>u]
area = sapply(ulocmax,auc,u=u,data=sdata); area = area[!is.na(area)]
#### Theoretical vs. Numerical values
vc = seq(0,30,0.5)
vGv = sapply(vc,Gv,lam1=lam1,lam2=lam2,u=u)
plot(vc,vGv,xlab="Excursion Area",ylab="Tail Probability",type="l",lwd=2,
     main=paste("Tail Probability of Excursion Area: (u = ",u,")"))
abline(h=0,lty=2,lwd=2,col="grey")
lines(x = ecdf.tail(area)[,1],y=ecdf.tail(area)[,2],col="blue",lwd=2)
lines(x = ecdf.tail(area_low)[,1],y=ecdf.tail(area_low)[,2],col="red",lwd=2)
legend("topright",legend=c("Theoretical value","Numerical value","Lower bound"),lwd=c(2,2,2),
       col=c("black","blue","red"),bty="n")
#### The distribution of Excursion area (Expectation form) is correct !!!
##############################################################
#### 2D
##############################################################
F2x = function(l,u,n=1000){
  # denominator
  F1 = function(x){
    Lam = matrix(rnorm(n*2),ncol=2)
    indicator = rep(NA,n)
    for(i in 1:n) indicator[i] = ifelse(Lam[i,2]>Lam[i,1] & Lam[i,2]<x/sqrt(2),1,0)
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2))*indicator)*dnorm(x)
  }
  # numerator
  F2 = function(x,l){
    Lam = matrix(rnorm(n*2),ncol=2)
    indicator1 = indicator2 = rep(NA,n)
    for(i in 1:n){
      indicator1[i] = ifelse(Lam[i,2] > Lam[i,1] & Lam[i,2] < x/sqrt(2),1,0)
      indicator2[i] = ifelse(4*pi^2*(x-u)^2 >= l^2*(Lam[i,2]-x/sqrt(2))*(Lam[i,1]-x/sqrt(2))/(2*4^4),1,0)
    } 
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2))*indicator1*indicator2)*dnorm(x)
  }
  vx = seq(u,20,0.2) 
  temp2 = trapz(vx,sapply(vx,F2,l=l))
  temp1 = trapz(vx,sapply(vx,F1))
  return(temp2/temp1)
}
####
L = 3000
h = seq(100,L,500)
sigma = seq(1,by=2,length.out=length(h))
signal = rep(0,L);region = vector()
for(i in 1:length(h)){
  c = 3*sqrt(2*pi)*sigma[i]
  temp = (h[i]-3*sigma[i]):(h[i]+3*sigma[i])
  signal[temp] = c*dnorm(temp,h[i],sigma[i])
  region = append(region,temp)
}
noise = rnorm(length(signal),0,1)
data = signal + noise
gamma = 4
sdata = smth.gau(data,gamma)
sdata = sdata*sqrt(sqrt(pi)*2*gamma) # standardize
ts.plot(sdata)
abline(h=u,lwd=2,lty=2,col="red")
abline(v=h,lwd=2,lty=2,col="blue")

u = 2.5
locmax = which.peaks(sdata)
ulocmax = locmax[sdata[locmax]>u]
width = sapply(lapply(ulocmax,search,u=u,sdata=sdata),sum)
tpid = which(ulocmax %in% region)

pval = sapply(width,Fx,lam1=lam1,lam2=lam2,u=u)
pthresh = fdrBH(pval,0.05)
pid = which(pval<=pthresh) 

Fdr(pid,tpid[-1],b=0)

###########################################################
u= 3
nr = 300; nc = 300   # 10000 * 10000
noise = matrix(rnorm(nr*nc),nr,nc)
sdata = smoothie::kernel2dsmooth(noise,"gauss",sigma=4,nx=nr,ny=nc)
sdata = sdata/0.0704  # 0.0704
# sigma = sqrt(1/(4*pi*gamma^2)) # gamma is the sigma in kernel2dsmooth
uloc = which(sdata>=u,arr.ind=TRUE)
plot(x=uloc[,2],y=uloc[,1],pch=20,xlab="column",ylab="row")
# 3D plot 
rgl::persp3d(x=1:nr,y=1:nc,z=sdata,phi=30,ltheta=120,col="lightblue",
        shade=0.75,xlab="X",ylab="Y",zlab="smoothed noise")
sqdf = data.frame(x=c(1,1,nr,nr,1),
                  y=c(1,nc,nc,1,1),
                  z=rep(u,5))
rgl::polygon3d(sqdf$x,sqdf$y,sqdf$z,coord=c(1,2),alpha=0.5,color="purple",add=T)
# heatmap 
sdata1 = sdata
sdata1[sdata1<u] = 0
fields::image.plot(x=1:nr,y=1:nc,sdata1,xlab="X",ylab="Y")
############
cross.area.cluster = function(data,u){
  x = data >= u
  nc = ncol(data); nr = nrow(data)
  find.line = function(x){
    # x is the data which contains only TRUE or FALSE
    columns = which(colSums(x) > 0) # non-False columns
    f = function(column){
      runs = rle(x[,column]) # consecutive repeat 
      runs_true = which(runs$value==TRUE)
      end = cumsum(runs$length)[runs_true]
      start = cumsum(runs$length)[runs_true-1] + 1
      # bottom element is TRUE
      if (0 %in% (runs_true-1)) start = c(1,start) 
      cbind(column,start,end)
    }
    as.data.frame(do.call("rbind",lapply(columns,f)))
  }
  ##
  find.cluster = function(x){
    # x is the output of function find.line
    match.line = function(x1,x2) {
      for(i in 1:nrow(x2)){
        ind = ifelse(!(x1$start > x2[i,]$end | x1$end < x2[i,]$start),i,0) 
        if(ind!=0) break 
      } 
      return(x2[ind,])
    }
    find.line.neighbor = function(y){
      # given a line, find all its neighboring lines
      cluster = y
      column = y$column + 1
      while(column %in% x$column){
        px = x[x$column==column,]
        y = match.line(y,px)
        if(nrow(y)==0) break # no match
        cluster = rbind(cluster,y)
        column = column + 1
      }
      return(cluster)
    }
    out = NULL; group = 0 
    columns = unique(x$column)
    for(j in columns){
      pat = x[x$column==j,]
      if(nrow(pat) == 0) next
      for(k in 1:nrow(pat)){
        group = group + 1
        out = rbind(out,cbind(find.line.neighbor(pat[k,]),group))
      }
      x = dplyr::anti_join(x,out[,1:3],by=colnames(x)) # delete matched lines
    }
    return(out)
  }
  ##
  line2point = function(x){
    # x is line-data 
    temp = rbind(cbind(x$column,x$end),cbind(rev(x$column),rev(x$start)))
    colnames(temp) = c("column","row") ; temp = as.data.frame(temp)
    return(cbind(temp[!duplicated(temp),],group=x$group[1]))
  }
  ##
  cross.area = function(x){
    # x is line-data 
    if("start" %in% colnames(x)) x = line2point(x)
    area = ifelse(nrow(x)<3,0,abs(polyarea(x[,1],x[,2])))
    return(data.frame(group=x$group[1],area=area))
  }
  ##
  outer = function(x){
    # x is line-data 
    t1 = x[1,]; t1$column = t1$column - 1
    t2 = x[nrow(x),]; t2$column = t2$column + 1
    x$start = x$start -1; x$end = x$end + 1
    temp = rbind(t1,x,t2)
    ## boundary issues
    temp = temp[!(temp$column %in% c(0,nc+1)),] #column boundary
    temp$start[temp$start==0] = 1; temp$start[temp$start==(nr+1)] = nr
    temp$end[temp$end==(nr+1)] = nr # row boundary
    return(temp)
  }
  ##
  fit.point = function(x){
    inner_point = line2point(x)
    outer_point = line2point(outer(x))
    euc.dist = function(x1, x2) sqrt(sum((x1 - x2)^2))
    interpolation = function(x){
      f1 = function(x) inner_point[apply(inner_point,1,euc.dist,x2=x) <= 1,]
      f2 = function(x1,y1,x2,y2) ((x2-x1)*u + x1*y2-x2*y1)/(y2-y1) # linear interpolation
      inner_match = f1(x)
      inpl_value = NULL; x = as.numeric(x)
      for(i in 1:nrow(inner_match)){
        t = inner_match[i,]; t = as.numeric(t) 
        iden_col = (x[1] == t[1])
        if(euc.dist(t,x)==0) inpl_value = rbind(inpl_value,x[1:2])
        else{
          temp = ifelse(iden_col,f2(x[2],data[x[2],x[1]],t[2],data[t[2],t[1]]),f2(x[1],data[x[2],x[1]],t[1],data[t[2],t[1]]))
          if(iden_col) inpl_value = rbind(inpl_value,cbind(x[1],temp))
          else inpl_value = rbind(inpl_value,cbind(temp,x[2]))
        }
      }
      colnames(inpl_value) = c("column","row")
      return(as.data.frame(inpl_value,stringsAsFactors = FALSE))
    } 
    if(nrow(x)==1) inner_point # single line cluster
    else as.data.frame(cbind(do.call(rbind,apply(outer_point,1,interpolation)),group=x$group[1]))
  }
  ## 
  clusters = find.cluster(find.line(x))
  outer_point=do.call("rbind",lapply(split(clusters,clusters$group),outer))
  fit_point=do.call("rbind",lapply(split(clusters,clusters$group),fit.point))
  area_lower = do.call("rbind",lapply(split(clusters,clusters$group),cross.area))
  area_upper = do.call("rbind",lapply(split(clusters,clusters$group),function(x) cross.area(outer(x))))
  area_fit = do.call("rbind",lapply(split(clusters,clusters$group),function(x) cross.area(fit.point(x))))
  return(list(cluster = clusters, outer_point = outer_point, fit_point = fit_point,
              area_lower = area_lower, area_upper = area_upper, area_fit = area_fit))
}
##
u = 3
xx = cross.area.cluster(sdata,u)
area.lower = xx$area_lower$area
area.upper = xx$area_upper$area 
area.fit = xx$area_fit$area

xlim = 0:50
plot(x = xlim,y=mean_tharea_u4,ylim=c(0,1),lwd=2,type="l",
     main = paste("u = ",u,seq=""),col = "black",xlab="area",ylab="")
lines(x = sort(area.lower),y=ecdf.tail(area.lower)[,2],lwd=2,col="blue",lty=2)
lines(x = sort(area.upper),y=ecdf.tail(area.upper)[,2],lwd=2,col="red",lty=2)
lines(x = sort(area.fit),y=ecdf.tail(area.fit)[,2],lwd=2,col="green",lty=1)
abline(h=0.1,lwd=2,col="grey")
legend("bottomleft",legend=c("theoretical","lower","upper","fitted"),
       lwd=rep(2,4),lty=c(1,2,2,1),col=c("black","blue","red","green"),bty="n")

plot(outer_fill1[,1:2],pch=25,bg="red")
points(inner_fill1[,1:2],pch=24,bg="black")
points(fit_fill1[,1:2],pch=21,bg="blue")
legend("topleft",legend = c("outer","inner","fitted"),pch=c(25,24,21),col=c("red","black","blue"),bty="n")
##
getVolume=function(df) {
  if(nrow(df) < 3 | length(unique(df$column))==1 | length(unique(df$row))==1) out = NA
  else{
    #find triangular tesselation of (x,y) grid
    res=geometry::delaunayn(as.matrix(df[,c("column","row")]),full=TRUE,options="Qz")
    #calulates sum of truncated prism volumes
    out = sum(mapply(function(triPoints,A) A/3*sum(df[triPoints,"z"]),
                     split.data.frame(res$tri,seq_along(res$areas)),
                     res$areas))
  }
  return(as.data.frame(cbind(group = df$group[1],volume = out,peakh = max(df$z))))
}
# sapply(split(x1,x1$group),getVolume)

## test example 
df1 <- data.frame(x=c(2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
                  y=c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3),
                  z=c(0,2,0,4,6,7,3,2,1,2,7,8,9,4,2))

rgl::persp3d(x=2:6,y=1:3,z=matrix(df1$z,nr=5,byrow=T),phi=30,ltheta=120,col="lightblue",
             shade=0.75,xlab="X",ylab="Y",zlab="smoothed noise")
sqdf = data.frame(x=c(1,1,nr,nr,1),
                  y=c(1,nc,nc,1,1),
                  z=rep(u,5))
rgl::polygon3d(sqdf$x,sqdf$y,sqdf$z,coord=c(1,2),alpha=0.5,color="purple",add=T)

# interpolation for the fitted points
# step1: fill all the points inside the boundary (fitted)
# step2: interpolation
# step3: calculate volume
# step4: calculate outer and inner volume
fill.point = function(x,data,isfit=F,clusters=NULL){
  # fill all the points that are inside the boundary
  # isfit: logic value indicates if x is fitted points
  if(isfit & is.null(clusters)) stop("clusters cannot be null for fited points") 
  if(isfit) {
     x0 = x
     x = clusters[clusters$group==x0$group[1],]}
  t = nrow(x); points = data.frame()
    for(i in 1:t){
      temp = x[i,]$start:x[i,]$end
      points = rbind(points,cbind(rep(x[i,]$column,length(temp)),temp))
    }
  colnames(points) = c("column","row")
  out = cbind(points,z = data[cbind(points$row,points$column)]-u,group=x$group[1])
  if(isfit) {
    out = rbind(out,cbind(column=x0$column,row=x0$row,z=0,group=x0$group[1]))
    out[order(out$column),]
  }
  return(out)
}
inner_fill=do.call("rbind",lapply(split(xx$cluster,xx$cluster$group),fill.point,data = sdata))
inner_vol = do.call("rbind",lapply(split(inner_fill,inner_fill$group),getVolume)); inner_vol = inner_vol[!is.na(inner_vol$volume),]
outer_fill=do.call("rbind",lapply(split(xx$outer_point,xx$outer_point$group),fill.point,data = sdata))
outer_vol = do.call("rbind",lapply(split(outer_fill,outer_fill$group),getVolume)); outer_vol = outer_vol[!is.na(outer_vol$volume),]
fit_fill = do.call("rbind",lapply(split(xx$fit_point,xx$fit_point$group),fill.point,data = sdata,isfit=T,clusters=xx$cluster))
fit_vol = do.call("rbind",lapply(split(fit_fill,fit_fill$group),getVolume)); fit_vol = fit_vol[!is.na(fit_vol$volume),]
#
xlim = 0:40
plot(x = xlim,y=mean_volume_u4[1:41],ylim=c(0,1),lwd=2,type="l",
     main = paste("u = ",u,seq=""),col = "black",xlab="volume",ylab="")
lines(x = sort(inner_vol),y=ecdf.tail(inner_vol)[,2],lwd=2,col="blue",lty=2)
lines(x = sort(outer_vol),y=ecdf.tail(outer_vol)[,2],lwd=2,col="red",lty=2)
lines(x = sort(fit_vol),y=ecdf.tail(fit_vol)[,2],lwd=2,col="green",lty=1)
abline(h=0.1,lwd=2,col="grey")
legend("topright",legend=c("theoretical","lower","upper","fitted"),
       lwd=rep(2,4),lty=c(1,2,2,1),col=c("black","blue","red","green"),bty="n")

#### explanation of fitted curve and theoretical values
library(raster)
heatp = function(x,xmn=0,ymn=0,xmx=1,ymx=1,nr=100,nc=100){
  x1 = as.data.frame(x)
  refgrid = raster(nrows=nr,ncols=nc,xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx)
  r = raster(refgrid); r[] = 0
  tab = table(cellFromXY(refgrid,x1));
  r[as.numeric(names(tab))] <- tab
  d <- data.frame(coordinates(r), count=r[])
  d = d[order(d$x,d$y),]
  #d[d$count>=50,"count"] = 50
  d$count = d$count/sum(d$count)
  fields::image.plot(x=unique(d$x),y=unique(d$y),matrix(d$count,nr,nc,byrow=T),
                     xlab="peak height",ylab="excursion",main="simulated data: u = 3")
}
#plot(inner_vol$peakh,inner_vol$volume)

################# 
## Joint distribution (peak height and excursion area)
## area
peak_area = data.frame(peakh = sapply(split(fit_fill,fit_fill$group),function(x) max(x$z)),
                       area = area.fit)
heatp(peak_area,xmx=2.2,ymx=155)
dev.off()
f1x = function(x) log(1-x)
f2x = function(x) -x
x = seq(-10,0.999,length.out=1000)
plot(x,f1x(x),type="l",lwd=2,col="red")
lines(x,f2x(x),col="blue",lwd=2)

joint.F2x = function(l,h,u,n=1000){
  # denominator
  F1 = function(x){
    Lam = matrix(rnorm(n*2),ncol=2)
    indicator = rep(NA,n)
    for(i in 1:n) indicator[i] = ifelse(Lam[i,2]>Lam[i,1] & Lam[i,2]<x/sqrt(2),1,0)
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2))*indicator)*dnorm(x)
  }
  # numerator
  F2 = function(x,l){
    Lam = matrix(rnorm(n*2),ncol=2)
    indicator1 = indicator2 = rep(NA,n)
    for(i in 1:n){
      indicator1[i] = ifelse(Lam[i,2] > Lam[i,1] & Lam[i,2] < x/sqrt(2),1,0)
      indicator2[i] = ifelse(4*pi^2*(x-u)^2 >= l^2*(Lam[i,2]-x/sqrt(2))*(Lam[i,1]-x/sqrt(2))/(2*4^4),1,0)
    } 
    mean((Lam[,2]-Lam[,1])*(Lam[,1]-x/sqrt(2))*(Lam[,2]-x/sqrt(2))*indicator1*indicator2)*dnorm(x)
  }
  vx = seq(u,20,0.2) 
  temp2 = trapz(vx+h,sapply(vx,F2,l=l))
  temp1 = trapz(vx,sapply(vx,F1))
  return(temp2/temp1)
}
#############
peak_width = data.frame(peakh = sdata[ulocmax]-u, width=width)
xmx = round(max(sdata[ulocmax]-u),1)
ymx = round(max(width),1)
heatp(peak_width,xmx=xmx,ymx=ymx)
# plot(peak_width)
x = seq(0,xmx,length.out=100)
y = seq(0,ymx,length.out=100)
fxh = function(h,x,lam1,lam2,u){
  if(lam1<=0 | lam2<=0) stop("lambda1 and lambda2 must be postive numbers")
  det = lam2-lam1^2
  den = lam2*x^4-16*lam1*x^2+64
  if(det<=0) stop("lambda2 is smaller than lambda1 square")
  uh = u + h
  num = lam2*x^2-8*lam1
  p1 = sqrt(lam2/(2*pi))*(1-pnorm(sqrt(lam2/det)*uh)) +
    lam1*dnorm(u)*(1-pnorm(-lam1*uh/sqrt(det)))
  p2 = 128*x*det/den^(3/2)*dnorm(8*u/sqrt(den)) # correct version
  p3 = (1+u^2*num^2/(den*det))*(1-pnorm(num*uh/(sqrt(det*den))))
  p4 = uh*num/(sqrt(det*den))*dnorm(uh*num/sqrt(det*den))
  return(1/p1*p2*(p3-p4))
}
z = outer(x,y,fxh,lam1,lam2,u)
fields::image.plot(x=x,y=y,z=z)




