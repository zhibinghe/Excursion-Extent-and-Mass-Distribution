####************* Required packages
library(dSTEM) # installed from github
library(pracma)
library(fields)

####************* Random seed
set.seed(2021)

#### Plotting settign
width = 8
height = 2

####************* Main functions

#### compute tail probability based on empirical cdf
ecdf.tail = function(x){
  fecdf = ecdf(x)
  x = sort(x)
  cbind(x,1 - fecdf(x))
}

#### search endpoints of f(x)=u
search = function(x,u,data){
  # x: location of a local maximum/minimum 
  # left-direction search
  lcount = 1
  while(data[x-lcount]>u & lcount < x-1){ lcount = lcount + 1}
  lcount = ifelse(lcount==x-1,NA,lcount)
  lcount = lcount - 1 + (data[x-lcount+1]-u)/(data[x-lcount+1]-data[x-lcount])
  # right-direction search
  rcount = 1
  while(data[x+rcount]>u & rcount < length(data)-x){ rcount = rcount + 1}  
  rcount = ifelse(rcount==length(data)-x,NA,rcount)
  rcount = rcount -1 + (data[x+rcount-1]-u)/(data[x+rcount-1]-data[x+rcount])
  return(c(lcount,rcount))
}

#### calculate excursion area for each local maximum/minimum
auc = function(x,u,data){
  # x: location of a local maximum/minimum
  endpoints = search(x,u,data)
  temp = (x-round(endpoints[1])):(x+round(endpoints[2]))
  traparea = trapz(x=temp,y=data[temp])- u*(tail(temp,1)-temp[1])
  traparea = ifelse(is.na(traparea),0,traparea)
}

#### 2D corss-surface for 'f(x,y) = u' estimate
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

#### fill all the points that are inside the boundary
fill.point = function(x,data,isfit=F,clusters=NULL){
  # x is boundary 
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
  if(!isfit) {
    out = rbind(out,cbind(column=x0$column,row=x0$row,z=0,group=x0$group[1]))
    out[order(out$column),]
  }
  return(out)
}

#### calculate excursion volume
getVolume=function(df) {
  # df is a dataframe from funtion fill.point
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
