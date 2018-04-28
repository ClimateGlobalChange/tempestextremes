jaccard_func<-function(x,y){
  xl<-length(x)
  yl<-length(y)
  sim<-length(which(x==y))
  j<-sim/(2*xl+2*yl-sim)
  print(j)
  return(j)
}
tl<-dim(pv_blob_sub)[3]
j=0
for (t in 1:tl){
  pv_vec<-as.vector(pv_blob_sub[,,t])
  gh_vec<-as.vector(gh_blob_sub[,,t])
  j<-j+jaccard_func(pv_vec,gh_vec)
}
j<-j/tl

#Random dist
a<-rnorm(1000,3,0.25)
b<-rnorm(1000,3.1,0.3)
c<-rnorm(1000,2.5,0.75)
d<-rnorm(1000,3,0.3)
plot(density(a),lwd=3,col="red",xlim=c(0,4))
lines(density(b),lwd=3,col="blue")
lines(density(c),lwd=3,col="green")
lines(density(d),lwd=3,col="purple")


