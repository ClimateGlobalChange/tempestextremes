#Create new vectors by rearranging the elements of Y1
# and Y2 and finding the difference in the means as before
# If Y1 and Y2 come from the same population, 
# then the original difference between the means should not be
# an extreme value in the distribution of differences
# If it does have an extreme value, then Y1 and Y2 are
# likely different populations
perm.diff<-function(Y1,Y2){
  Y<-c(Y1,Y2)
  l2<-length(Y)/2
  Ysamp<-sample(Y,length(Y))
  Y1p<-Ysamp[1:l2]
  Y2p<-Ysamp[(l2+1):length(Y)]
  dp<-mean(Y1p)-mean(Y2p)
  return(dp)
}

calc.perm.p<-function(Y1,Y2,n){
  Y1s<-sample(Y1,n)
  Y2s<-sample(Y2,n)
  d1<-mean(Y1s)-mean(Y2s)
  
  U1<-replicate(999,perm.diff(Y1s,Y2s))
  diff1<-sort(c(U1,d1))
  # how many values of U are less than original difference?
  # If this value is very small or large, then d is an extreme
  # Therefore, we can reject the null that the two populations
  # are from the same distribution
  U1.low<-U1[U1<d1]

  #Count the number in the tail, i.e. the number of
  #values that are further from the median than d
  n1.tail<-ifelse(length(U1.low)<=length(diff1)/2,length(U1.low),length(diff1)-length(U1.low))
 #print(sprintf("length of tail: %f length of vector: %f",n1.tail,length(diff1)))
  #Double the number to account for both tails
  p1<-2*n1.tail/length(diff1)
  return(p1)
}

calc.sample.size<-function(arr1,arr2,lons,lats){
  coords<-expand.grid(lons,lats)
  tsize<-dim(arr1)[3]
  lonsize<-dim(arr1)[2]
  sampsize<-c()
  for (t in 1:tsize){
    vec1<-as.vector(arr1[,,t])
    vec2<-as.vector(arr2[,,t])
    #Express the data as a surface
    surf1<-surf.ls(0,coords[,1],coords[,2],vec1)
    surf2<-surf.ls(0,coords[,1],coords[,2],vec2)
    r1<-correlogram(surf1,lonsize,plotit=FALSE)
    r2<-correlogram(surf2,lonsize,plotit=FALSE)
    
    nr1r2<-r1$cnt*r1$y*r2$y
    if (!is.na(sum(nr1r2))){
      sr.hat<-sum(nr1r2)/length(nr1r2)^4
      ne.hat<-1+1/sr.hat
      sampsize<-c(sampsize,ne.hat)
    }
  }
  return(mean(na.omit(sampsize)))
}


W.calc<-function(mat){
  n.tot<-sum(mat)
  mat<-mat+1
  theta.hat<-(mat[1,1]*mat[2,2])/(mat[1,2]*mat[2,1])
  sig2.hat<-sum(1/(mat/n.tot))
  W<-n.tot*log(theta.hat)^2/sig2.hat
  return(W)
}

cont.table<-function(vec1,vec2){
  n.mat<-matrix(0,nrow=2,ncol=2)
  colnames(n.mat)<-c("V1","V2")
  rownames(n.mat)<-c("Pres","Abs")
  
  n.mat[1,1]<-length(which(vec1>0))
  n.mat[1,2]<-length(which(vec2>0))
  n.mat[2,1]<-length(vec1)-n.mat[1,1]
  n.mat[2,2]<-length(vec2)-n.mat[1,2]
  return(n.mat)
}

pres.abs.table<-function(vec1,vec2){
  presabs<-ifelse(vec1==0,ifelse(vec2==0,"None","V2"),ifelse(vec2==0,"V1","Both"))
  n.mat<-matrix(0,nrow=2,ncol=2)
  colnames(n.mat)<-c("Abs","Pres")
  rownames(n.mat)<-c("Abs","Pres")
  
  n.mat[1,1]<-length(which(presabs=="None"))
  n.mat[1,2]<-length(which(presabs=="V1"))
  n.mat[2,1]<-length(which(presabs=="V2"))
  n.mat[2,2]<-length(which(presabs=="Both"))
  return(n.mat)
}

W.adj<-function(pv.vec,z.vec,lon_seq,lat_seq){
  n.mat<-matrix(0,nrow=2,ncol=2)
  colnames(n.mat)<-c("V1","V2")
  rownames(n.mat)<-c("Pres","Abs")
  
  n.mat[1,1]<-length(which(pv.vec>0))
  n.mat[1,2]<-length(which(z.vec>0))
  n.mat[2,1]<-length(pv.vec)-n.mat[1,1]
  n.mat[2,2]<-length(z.vec)-n.mat[1,2]
  
  chi1<-chisq.test(n.mat,simulate.p.value = TRUE)
  
  n.tot<-sum(n.mat)
  p1<-(n.mat[1,1]+n.mat[1,2])/n.tot
  p2<-(n.mat[2,1]+n.mat[2,2])/n.tot
  q1<-(n.mat[1,1]+n.mat[2,1])/n.tot
  q2<-(n.mat[2,1]+n.mat[2,2])/n.tot
  
  
  coords<-expand.grid(lon_seq,lat_seq)
  
  pv.surf<-surf.ls(0,coords[,1],coords[,2],pv.vec)
  pv.vgr<-variogram(pv.surf,length(lon_seq),plotit = FALSE)
  CX.hat<-max(pv.vgr$y[1:50])-pv.vgr$y[1:50]
  nkX<-pv.vgr$cnt[1:50]
  
  z.surf<-surf.ls(0,coords[,1],coords[,2],z.vec)
  z.vgr<-variogram(z.surf,length(lon_seq),plotit = FALSE)
  CY.hat<-max(z.vgr$y[1:50])-z.vgr$y[1:50]
  nkY<-pv.vgr$cnt[1:50]
  
  lambda.sum <- sum(CX.hat * CY.hat * nkX)
  
  lambda <-(2 * n.tot* lambda.sum)/(n.mat[1,1]*n.mat[2,2])
  W<-W.calc(n.mat)
  W.adj <- chi1$statistic / (1 + lambda)
  p <- 1 - pchisq(W.adj, df = 1)
  return(c(W.adj,p,lambda))
}