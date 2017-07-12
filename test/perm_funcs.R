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
  n_sample<-sample(1:length(Y1),n)
  Y1s<-Y1[n_sample]
  Y2s<-Y2[n_sample]
  d1<-mean(Y1s)-mean(Y2s)
  #print(sprintf("The mean is %f",d1))
  U1<-replicate(9999,perm.diff(Y1s,Y2s))
  diff1<-sort(c(U1,d1))
  # how many values of U are less than original difference?
  # If this value is very small or large, then d is an extreme
  # Therefore, we can reject the null that the two populations
  # are from the same distribution
  U1.low<-U1[U1<d1]
  #print("Vector less than mean:")
  #print(U1.low)
  #print(sprintf("The length of this vector is %d",length(U1.low)))
  #Count the number in the tail, i.e. the number of
  #values that are further from the median than d
  n1.tail<-ifelse(length(U1.low)<=length(diff1)/2,length(U1.low),length(diff1)-length(U1.low))
  #print(sprintf("length of tail: %f length of vector: %f",n1.tail,length(diff1)))
  #Double the number to account for both tails
  p1<-2*n1.tail/length(diff1)
  return(p1)
}


block.perm.calc<-function(arr1,arr2,block_edge,lon_seq,lat_seq){
  
  nlonbox<-as.integer(length(lon_seq)/block_edge)
  nbx<-nlonbox*block_edge
  rbx<-length(lon_seq)-nbx
  nlatbox<-as.integer(length(lat_seq)/block_edge)
  nblocks<-nlonbox*nlatbox
  nby<-nlatbox*block_edge
  rby<-length(lat_seq)-nby
  
  #trim the edges
  nsubleft<-0
  nsubright<-0
  nsubbot<-0
  nsubtop<-0
  if (rbx>0){
    nsubleft<-as.integer(rbx/2)
    nsubright<-rbx-nsubleft
  }
  if (rby>0){
    nsubbot<-as.integer(rby/2)
    nsubtop<-rby-nsubbot
  }
  loninds<-seq((1+nsubleft),(length(lon_seq)-nsubright))
  latinds<-seq((1+nsubbot),(length(lat_seq)-nsubtop))
  new1<-arr1[loninds,latinds]
  new2<-arr2[loninds,latinds]
  
  B1<-array(NA,c(nblocks,block_edge,block_edge))
  B2<-array(NA,c(nblocks,block_edge,block_edge))
  
  
  b<-1
  for (m in 1:nlonbox){
    for (n in 1:nlatbox){
      l1<-seq(m-1,m)*block_edge
      l1[1]<-l1[1]+1
      l2<-seq(n-1,n)*block_edge
      l2[1]<-l2[1]+1
      B1[b,,]<-new1[l1[1]:l1[2],l2[1]:l2[2]]
      B2[b,,]<-new2[l1[1]:l1[2],l2[1]:l2[2]]
      b<-b+1
    }
  } 
  #Resampled blocks
  i.samp<-sample(1:nblocks,replace=TRUE)
  C1<-array(NA,c(nbx,nby))
  C2<-array(NA,c(nbx,nby))
  x<-1
  for (m in 1:nlonbox){
    for (n in 1:nlatbox){
      l1<-seq(m-1,m)*block_edge
      l1[1]<-l1[1]+1
      l2<-seq(n-1,n)*block_edge
      l2[1]<-l2[1]+1
      C1[l1[1]:l1[2],l2[1]:l2[2]]<-B1[i.samp[x],,]
      C2[l1[1]:l1[2],l2[1]:l2[2]]<-B2[i.samp[x],,]
      x<-x+1
    }
  }
  
  d<-mean(C1)-mean(C2)
  return(d)
}

calc.perm.block<-function(arr1,arr2,block_edge,lon_seq,lat_seq){
  d1<-mean(arr1)-mean(arr2)
  drep<-replicate(999,block.perm.calc(arr1,arr2,block_edge,lon_seq,lat_seq))
  diff1<-sort(c(drep,d1))
  d.low<-drep[drep<d1]
  n1.tail<-ifelse(length(d.low)<=length(diff1)/2,length(d.low),length(diff1)-length(d.low))
  #print(sprintf("length of tail: %f length of vector: %f",n1.tail,length(diff1)))
  #Double the number to account for both tails
  p1<-2*n1.tail/length(diff1)
  return(p1)
}

# calc.perm.p<-function(Y1,Y2,n){
#   Y1s<-sample(Y1,n)
#   Y2s<-sample(Y2,n)
#   d1<-mean(Y1s)-mean(Y2s)
#   print(sprintf("The mean is %f",d1))
#   U1<-replicate(9999,perm.diff(Y1s,Y2s))
#   diff1<-sort(c(U1,d1))
#   # how many values of U are less than original difference?
#   # If this value is very small or large, then d is an extreme
#   # Therefore, we can reject the null that the two populations
#   # are from the same distribution
#   U1.low<-U1[U1<d1]
#   print("Vector less than mean:")
#   print(U1.low)
#   print(sprintf("The length of this vector is %d",length(U1.low)))
#   #Count the number in the tail, i.e. the number of
#   #values that are further from the median than d
#   n1.tail<-ifelse(length(U1.low)<=length(diff1)/2,length(U1.low),length(diff1)-length(U1.low))
#  print(sprintf("length of tail: %f length of vector: %f",n1.tail,length(diff1)))
#   #Double the number to account for both tails
#   p1<-2*n1.tail/length(diff1)
#   return(p1)
# }

#Effective sample size
calc.sample.size<-function(vec1,vec2,lons,lats){
  coords<-expand.grid(lons,lats)
  lonsize<-length(lons)
  surf1<-surf.ls(0,coords[,1],coords[,2],vec1)
  surf2<-surf.ls(0,coords[,1],coords[,2],vec2)
  r1<-correlogram(surf1,lonsize,plotit=FALSE)
  r2<-correlogram(surf2,lonsize,plotit=FALSE)
  
  nr1r2<-r1$cnt*r1$y*r2$y
  sr.hat<-sum(nr1r2)/length(nr1r2)^4
  ne.hat<-1+1/sr.hat
  return(as.integer(ne.hat))
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