load("~/block_r_data/eur_cold_spell.RData")
#load("~/block_r_data/stats_cold_spell.RData")
library(ade4)
library(spdep)

#Comparing distances using the mantel test
#1987-01-02 0Z has overlaps between all 3 + an overlap of PV and Z
pv_sub<-pv_blob_sub[,,128]
pv_nonz<-which(pv_sub>0)
pv_sub[pv_nonz]<-pv_sub[pv_nonz]/pv_sub[pv_nonz]

z_sub<-z_blob_sub[,,128]
z_nonz<-which(z_sub>0)
z_sub[z_nonz]<-z_sub[z_nonz]/z_sub[z_nonz]

gh_sub<-gh_blob_sub[,,128]
gh_nonz<-which(gh_sub>0)
gh_sub[gh_nonz]<-gh_sub[gh_nonz]/gh_sub[gh_nonz]

#Create a sample of the two matrices
set.seed(1)
#total number of indices
ni<-length(pv_sub)
ind<-sample(seq(1,ni),2000,replace=FALSE)

pv_sample<-pv_sub[ind]
z_sample<-z_sub[ind]
g_sample<-gh_sub[ind]
#first cor
c1<-cor.test(pv_sample,z_sample,method="kendall")

#Now scramble Z
zval<-c(c1$statistic)
pval<-c(c1$p.value)
tval<-c(c1$estimate)
for (x in 1:1000){
  z_sample2<-sample(z_sample,length(z_sample),replace=FALSE)
  c2<-cor.test(pv_sample,z_sample2,method="kendall")
  zval<-c(zval,c2$statistic)
  pval<-c(pval,c2$p.value)
  tval<-c(tval,c2$estimate)
}

#first cor
c3<-cor.test(g_sample,z_sample,method="kendall")

#Now scramble Z
zval2<-c(c3$statistic)
pval2<-c(c3$p.value)
tval2<-c(c3$estimate)
for (x in 1:1000){
  z_sample2<-sample(z_sample,length(z_sample),replace=FALSE)
  c4<-cor.test(g_sample,z_sample2,method="kendall")
  zval2<-c(zval2,c4$statistic)
  pval2<-c(pval2,c4$p.value)
  tval2<-c(tval2,c4$estimate)
}
