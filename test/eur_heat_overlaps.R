load("~/block_r_data/eur_cold_spell.RData")
library("maps")
library("maptools")
library("mapproj")
#test overlaps with heat wave data

#All Z* has value of 100
nonz<-which(z_blob_sub>0)
z_blob2<-z_blob_sub
z_blob2[nonz]<-z_blob2[nonz]/z_blob2[nonz]
z_blob2<-z_blob2*100
#All PV* has value of 10
nonz<-which(pv_blob_sub>0)
pv_blob2<-pv_blob_sub
pv_blob2[nonz]<-pv_blob2[nonz]/pv_blob2[nonz]
pv_blob2<-pv_blob2*1000

#All ZG has value of 100
nonz<-which(gh_blob_sub>0)
gh_blob2<-gh_blob_sub
gh_blob2[nonz]<-gh_blob2[nonz]/gh_blob2[nonz]
gh_blob2<-gh_blob2*10000

overlaps_sub<-z_blob2+pv_blob2+gh_blob2
overlaps_nonz<-nonz<-overlaps_sub[which(overlaps_sub>0)]
#Counts of overlap/nonoverlap
overlaps_tab<-table(overlaps_nonz)
pct_overlaps<-overlaps_tab/sum(overlaps_tab)*100
pct_fmt<-sprintf("%2.2f",pct_overlaps)
labs<-c("Z*","PV*","PV*/Z*","ZG",
        "Z*/ZG","PV*/ZG","all")
labs_pct<-paste(labs,pct_fmt)
pie(overlaps_tab,
    labels=labs_pct,radius=1,
    col=c("blue","darkgreen","gold","purple",
          "magenta","darkorange","red"),
    main="Relative percentages of blocked area")

#Overlap of PV and GH
overlaps_PZA<-ifelse(overlaps_sub==1100,1,0)
#Overlap of Z and ZG
overlaps_ZZG<-ifelse(overlaps_sub==10100,1,0)
#Overlap of PV and ZG
overlaps_PZG<-ifelse(overlaps_sub==11000,1,0)
#Overlap of all 3
overlaps_all<-ifelse(overlaps_sub==11100,1,0)
hgt.cols<-colorRampPalette(c("purple","blue","cyan4","green","yellow", "orange","red","darkred"))(34)
for (t in first_ind:length(time_format)){
  #for(t in first_ind:first_ind+3){
  fname<-sprintf("~/Dropbox/overlaps/heat/heat_overlaps_%s_%02dZ_noT.png",time_format[t],time_hours[t])
  png(fname,height=600,width=800)
  
  map('world',xlim=c(-110,50),ylim=c(25,75),fill=TRUE)
  title(sprintf("Z500 %s %02dZ, Overlaps for PV*/Z* (yellow), PV*/ZG (orange), Z*/ZG (pink), all (red)",time_format[t],time_hours[t]))
  map.axes()
  contour(lon_seq,lat_seq,z_hgt_sub[,,t],levels=seq(4500,6100,50),add=TRUE,col=hgt.cols,lwd=2)
  contour(lon_seq,lat_seq,overlaps_PZA[,,t],levels=c(0,1),add=TRUE,col="gold",
          drawlabels=FALSE,lwd=4)
  contour(lon_seq,lat_seq,overlaps_PZG[,,t],levels=c(0,1),add=TRUE,col="darkorange",
          drawlabels=FALSE,lwd=4)
  contour(lon_seq,lat_seq,overlaps_ZZG[,,t],levels=c(0,1),add=TRUE,col="magenta",
          drawlabels=FALSE,lwd=4)
  contour(lon_seq,lat_seq,overlaps_all[,,t],levels=c(0,1),add=TRUE,col="red",
          drawlabels=FALSE,lwd=4)
  dev.off()
}