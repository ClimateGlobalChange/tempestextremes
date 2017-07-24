source("~/tempestextremes/test/sector_funcs.R")

rv_calc<-function(umat,vmat,dphi,dlambda,cosphi,ntime,nlat,nlon){
  invdphi<-1/(2*dphi)
  invdl<-1/(2*dlambda)
  dUdphi<-array(0,c(nlon,nlat,ntime))
  dVdl<-array(0,c(nlon,nlat,ntime))
  rv<-array(0,c(nlon,nlat,ntime))
  #partial derivatives!
  #U wrt phi
  for (b in 1:nlon){
    for (t in 1:ntime){
      dUdphi[b,1,t]<-(-umat[b,3,t]*cosphi[3]+4*umat[b,2,t]*cosphi[2]-3*umat[b,1,t]*cosphi[1])*invdphi
      dUdphi[b,nlat,t]<-(3*umat[b,nlat,t]*cosphi[nlat]-4*umat[b,nlat-1,t]*cosphi[nlat-1]+umat[b,nlat-2,t]*cosphi[nlat-2])*invdphi
      for (a in 2:(nlat-1)){
        dUdphi[b,a,t]<-(umat[b,a+1,t]*cosphi[a+1]-umat[b,a-1,t]*cosphi[a-1])*invdphi
      }
    }
  }
  #V wrt lambda
  for (a in 1:nlat){
    for (t in 1:ntime){
      dVdl[1,a,t]<-(vmat[2,a,t]-vmat[nlon,a,t])*invdl
      dVdl[nlon,a,t]<-(vmat[1,a,t]-vmat[nlon-1,a,t])*invdl
      for (b in 2:(nlon-1)){
        dVdl[b,a,t]<-(vmat[b+1,a,t]-vmat[b-1,a,t])*invdl
      }
    }
  }
  
  for (b in 1:nlon){
    for (a in 1:nlat){
      coef<-1/(6371000*cosphi[a])
      for (t in 1:ntime){
        rv[b,a,t]<-coef*(dVdl[b,a,t]-dUdphi[b,a,t])
      }
    }
  }
  return(rv)
}

fname="/Volumes/ExFAT_drive/ERA_files/ERA_2003/ERA_2003_08_vars.nc"
f<-open.nc(fname)
LEFT_BOUND=c(250, 30, 130, 290, 20, 120)
RIGHT_BOUND=c(50, 150, 270, 40, 140, 310)
MIN_LAT=c(25, 25, 25, -75, -75, -75)
MAX_LAT=c(75, 75, 75, -25, -25, -25)

lats<-read.nc(f)$lat
dphi<-abs(lats[1]-lats[2])
lons<-read.nc(f)$lon
dlambda<-abs(lons[1]-lons[2])
levs<-read.nc(f)$lev
timevar<-read.nc(f)$time


#Get array indices for bounds
li<-get_index(lons,LEFT_BOUND[1])
ri<-get_index(lons,RIGHT_BOUND[1])
bi<-get_index(lats,MIN_LAT[1])
ti<-get_index(lats,MAX_LAT[1])
lon_length<-ri+(length(lons)-li+1)
lat_length<-(bi-ti)+1
time_length<-length(timevar)
# 
lon_seq<-seq(-110,50)
lat_seq<-seq(25,75)
cosphi<-cos(lat_seq*pi/180)

u500<-array(NA,c(lon_length,lat_length,time_length))
u500[1:(lon_length-ri),,]<-read.nc(f)$U[li:length(lons),bi:ti,10,]
u500[(lon_length-ri+1):lon_length,,]<-read.nc(f)$U[1:ri,bi:ti,10,]
v500<-array(NA,c(lon_length,lat_length,time_length))
v500[1:(lon_length-ri),,]<-read.nc(f)$V[li:length(lons),bi:ti,10,]
v500[(lon_length-ri+1):lon_length,,]<-read.nc(f)$V[1:ri,bi:ti,10,]

close.nc(f)

z_hgt_sub<-array(NA,c(lon_length,lat_length,time_length))
fname="/Volumes/ExFAT_drive/ERA_files/ERA_2003/ERA_2003_08_vars_z500.nc"
f<-open.nc(fname)
z_hgt_sub[1:(lon_length-ri),,]<-read.nc(f)$Z[li:length(lons),bi:ti,]
z_hgt_sub[(lon_length-ri+1):lon_length,,]<-read.nc(f)$Z[1:ri,bi:ti,]
close.nc(f)

f<-open.nc("/Volumes/ExFAT_drive/ERA_files/ERA_blobs/ERA_2003_JJA_NA_blobs.nc")
timevar2<-read.nc(f)$time
tl<-length(timevar2)

time_format2<-as.Date(timevar2/24, origin="1800-01-01")
time_hours2<-timevar%%24




pv_blob_sub<-array(NA,c(lon_length,lat_length,tl))
pv_blob_sub[1:(lon_length-ri),,]<-read.nc(f)$PV_BLOB[li:length(lons),bi:ti,]
pv_blob_sub[(lon_length-ri+1):lon_length,,]<-read.nc(f)$PV_BLOB[1:ri,bi:ti,]
close.nc(f)
# 
# 
z_blob_sub<-array(NA,c(lon_length,lat_length,tl))
f<-open.nc("/Volumes/ExFAT_drive/ERA_files/ERA_blobs/ERA_2003_JJA_NA_Zblobs.nc")
z_blob_sub[1:(lon_length-ri),,]<-read.nc(f)$Z_BLOB[li:length(lons),bi:ti,]
z_blob_sub[(lon_length-ri+1):lon_length,,]<-read.nc(f)$Z_BLOB[1:ri,bi:ti,]
close.nc(f)
# 
gh_blob_sub<-array(NA,c(lon_length,lat_length,tl))
f<-open.nc("/Volumes/ExFAT_drive/ERA_files/ERA_blobs/ERA_2003_JJA_NA_GHGblobs.nc")
gh_blob_sub[1:(lon_length-ri),,]<-read.nc(f)$GHG_BLOB[li:length(lons),bi:ti,]
gh_blob_sub[(lon_length-ri+1):lon_length,,]<-read.nc(f)$GHG_BLOB[1:ri,bi:ti,]
close.nc(f)

rv<-rv_calc(u500,v500,dphi,dlambda,cosphi,time_length,lat_length,lon_length)
time_format<-as.Date(timevar/24, origin="1800-01-01")
time_hours<-timevar%%24

target_date<-"2003-08-07"

t<-which(time_format==target_date)
t<-t+1
t2<-which(time_format2==target_date)
t2<-t2+1
#
rv.cols<-colorRampPalette(c("blue","white","red"))(17)
rv.cols[8]<-rv.cols[9]

#for(t in first_ind:first_ind+3){
#fname<-sprintf("~/Dropbox/eur_heat_wave/tcol/eur_heat_%s_%02dZ.png",time_format[t],time_hours[t])
#png(fname,height=600,width=800)
filled.contour(lon_seq,lat_seq,rv[,,t],main=sprintf("Relative vorticity %s %02dZ",time_format[t],time_hours[t]),
               levels=seq(-4e-6,4e-6,0.5e-6),col=rv.cols,plot.axes={axis(1);axis(2); 
                 map('world',add=TRUE,xlim=c(-110,50),ylim=c(25,75));  
                 contour(lon_seq,lat_seq,z_hgt_sub[,,t],levels=seq(4500,6100,50),drawlabels=FALSE,add=TRUE);
                 contour(lon_seq,lat_seq,pv_blob_sub[,,t2],levels=c(0,1),add=TRUE,col="darkgreen",drawlabels=FALSE,lwd=4)
                 contour(lon_seq,lat_seq,z_blob_sub[,,t2],levels=c(0,1),add=TRUE,col="blue",drawlabels=FALSE,lwd=4)
                 contour(lon_seq,lat_seq,gh_blob_sub[,,t2],levels=c(0,1),add=TRUE,col="purple",drawlabels=FALSE,lwd=4)
                 #contour(lon_seq,lat_seq,rv[,,t],levels=seq(-4e-6,4e-6,0.5e-6),add=TRUE,drawlabels=FALSE,lwd=1)
                 })  
#dev.off()
