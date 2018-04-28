# library(RNetCDF)
# 
# #load an averaged data set for PV
# f<-open.nc("/Volumes/ExFAT_drive/ERA_files/ERA_avg/ERA_1980_2005_daily_avg.nc")
# lats<-read.nc(f)$lat
# #lats[46]=45
# lons<-read.nc(f)$lon
# #lons[226]=225
# loni<-75
# lati<-20
# day_in_year<-read.nc(f)$time
# #to get info: print.nc(f)
# aipv<-read.nc(f)$AIPV
# close.nc(f)
# aipv_sub<-aipv[loni,lati,]
# 
# 
# #data goes from 1980-2005
# #Loop through
# y<-seq(1980,2005)
# m<-seq(1,12)
# t<-c()
# ipv<-c()
# 
# for (a in y){
#   subdir<-sprintf("/Volumes/ExFAT_drive/ERA_files/ERA_%d/",a)
#   for (b in m){
#     fname<-sprintf("%sERA_%d_%02d_vars_integ.nc",subdir,a,b)
#     print(fname)
#     f<-open.nc(fname)
#     ft<-read.nc(f)$time
#     nt<-length(ft)
#     fipv<-read.nc(f)$IPV[loni,lati,]
#     if (b==2 & nt>28*4){
#       ft<-ft[1:(nt-4)]
#       fipv<-fipv[1:(nt-4)]
#     }
#     t<-c(t,ft)
#     ipv<-c(ipv,fipv)
#     close.nc(f)
#   }
# }
# 
# time_format<-as.Date(t/24, origin="1800-01-01")
# day_num<-as.numeric(strftime(time_format,format="%j"))
# 
# #Indices where new year starts
# ny<-length(ipv)/(365*4)
# nums_start<-seq(0,ny-1)*365*4+1
# nums_end<-seq(1,ny)*365*4
# day_num_3<-c()
# for (d in 1:365){
#   day_num_3<-c(day_num_3,rep(d,4))
# }
load("~/block_r_data/harmonics_data_1980-2005.RData")
plot(day_in_year,aipv_sub,type="l",ylim=c(0,8e-6),col="red",lwd=2,xlim=c(1,365),
     xlab="Day in year",ylab="Integrated PV",main="PV, averaged (black) with LTDM (blue), \n constant threshold (green) and 1 sd (purple)")

ipv_yavg<-rep(0,length(day_num_3))
avg_each<-rep(aipv_sub,each=4)
thresh<-rep(0,length(day_num_3))

for (i in 1:length(nums_start)){
  ipv_sub<-ipv[nums_start[i]:nums_end[i]]
  ipv_yavg<-ipv_yavg+ipv_sub
  idev<-ipv_sub-avg_each
  thresh<-thresh+(idev*idev)
  lines(day_num_3,ipv_sub,col="grey")
}
ipv_yavg<-ipv_yavg/ny
sd<-sqrt(thresh/length(nums_start))

daily_ipv<-rep(0,365)
sd_daily_pv<-rep(0,365)
for (i in 1:365){
  sub<-ipv_yavg[((i-1)*4+1):((i-1)*4+4)]
  daily_ipv[i]<-mean(sub)
  sub<-sd[((i-1)*4+1):((i-1)*4+4)]
  sd_daily_pv[i]<-mean(sub)
}

fft_sd_pv<-fft(sd_daily_pv)
fft_sd_pv[7:(length(fft_sd_pv)-6)]<-0
invfft_sd_pv<-fft(fft_sd_pv,inverse=TRUE)/length(fft_sd_pv)

fft_ipv<-fft(daily_ipv)
fft_ipv[7:(length(fft_ipv)-6)]<-0
invfft<-fft(fft_ipv,inverse=TRUE)/length(fft_ipv)

lines(day_in_year,invfft,lwd=2,col="blue")
lines(day_in_year,daily_ipv)
#lines(day_in_year,aipv_sub,col="red",lwd=2)
lines(day_in_year, invfft-invfft_sd_pv,col="purple",lwd=2)
lines(day_in_year, invfft-1.2e-6,col="darkgreen",lwd=2)

fft_ipv2<-fft(daily_ipv)
fft_ipv2[21:(length(fft_ipv2)-20)]<-0
invfft2<-fft(fft_ipv2,inverse=TRUE)/length(fft_ipv2)

#lines(day_in_year,invfft2,lwd=2,col="purple")

abline(v=59,lwd=2)
abline(v=151,lwd=2)
abline(v=243,lwd=2)
abline(v=334,lwd=2)

text((59+151)/2,6e-6,"MAM")
text((243+151)/2,6e-6,"JJA")
text((243+334)/2,6e-6,"SON")
text(59/2,6e-6,"DJF")
text((334+366)/2,6e-6,"DJF")

plot(day_in_year,sd_daily_pv,main="PV, standard deviation threshold (daily=black, FFT=blue)\n with constant threshold (green)",
     type="l",xlab="Day in year",ylab="PV (PVU)")
lines(day_in_year,invfft_sd_pv,col="blue",lwd=2)
abline(h=1.2e-6,lwd=2,col="darkgreen")

abline(v=59,lwd=2)
abline(v=151,lwd=2)
abline(v=243,lwd=2)
abline(v=334,lwd=2)

text((59+151)/2,180,"MAM")
text((243+151)/2,180,"JJA")
text((243+334)/2,180,"SON")
text(59/2,180,"DJF")
text((334+366)/2,180,"DJF")


###### Repeat for Z because I'm lazy

# ##load an averaged data set for PV
# f<-open.nc("/Volumes/ExFAT_drive/ERA_files/ERA_avg/ERA_1980_2005_Z500_avg.nc")
# lats<-read.nc(f)$lat
# #lats[46]=45
# lons<-read.nc(f)$lon
# #lons[226]=225
# day_in_year<-read.nc(f)$time
# #to get info: print.nc(f)
# az<-read.nc(f)$AVGZ
# close.nc(f)
# az_sub<-az[loni,lati,]
# 
# 
# #data goes from 1980-2005
# #Loop through
# y<-seq(1980,2005)
# m<-seq(1,12)
# t<-c()
# z<-c()
# 
# for (a in y){
#   subdir<-sprintf("/Volumes/ExFAT_drive/ERA_files/ERA_%d/",a)
#   for (b in m){
#     fname<-sprintf("%sERA_%d_%02d_vars_z500.nc",subdir,a,b)
#     print(fname)
#     f<-open.nc(fname)
#     ft<-read.nc(f)$time
#     nt<-length(ft)
#     fz<-read.nc(f)$Z[loni,lati,]
#     if (b==2 & nt>28*4){
#       ft<-ft[1:(nt-4)]
#       fz<-fz[1:(nt-4)]
#     }
#     t<-c(t,ft)
#     z<-c(z,fz)
#     close.nc(f)
#   }
# }
# 
# time_format<-as.Date(t/24, origin="1800-01-01")
# day_num<-as.numeric(strftime(time_format,format="%j"))
# 
# 
# fft_z<-fft(z)
# fft_z[20:length(fft_z)]<-0
# invfftz<-fft(fft_z,inverse=TRUE)/length(fft_z)
# 

#Indices where new year starts
ny<-length(z)/(365*4)
nums_start<-seq(0,ny-1)*365*4+1
nums_end<-seq(1,ny)*365*4

plot(day_in_year,az_sub,type="l",ylim=c(4800,5800),col="red",lwd=2,xlim=c(1,365),
     xlab="Day in year",ylab="Z",main="Z, averaged (black) with LTDM (blue), \n constant threshold (green) and 1 sd (purple)")

z_yavg<-rep(0,length(day_num_3))
avg_each<-rep(az_sub,each=4)
thresh<-rep(0,length(day_num_3))
for (i in 1:length(nums_start)){
  z_sub<-z[nums_start[i]:nums_end[i]]
  z_yavg<-z_yavg+z_sub
  iz<-z_sub-avg_each
  thresh<-thresh+(iz*iz)
  lines(day_num_3,z_sub,col="grey")
}
z_yavg<-z_yavg/ny
sd<-sqrt(thresh/length(nums_start))
#lines(day_num_3,z_yavg)

daily_z<-rep(0,365)
sd_daily_z<-rep(0,365)
for (i in 1:365){
  sub<-z_yavg[((i-1)*4+1):((i-1)*4+4)]
  daily_z[i]<-mean(sub)
  sub<-sd[((i-1)*4+1):((i-1)*4+4)]
  sd_daily_z[i]<-mean(sub)
}

fft_sd_z<-fft(sd_daily_z)
fft_sd_z[7:(length(fft_sd_z)-6)]<-0
invfft_sd_z<-fft(fft_sd_z,inverse=TRUE)/length(fft_sd_z)

#lines(day_in_year,az_sub,col="red",lwd=2)

fft_z<-fft(daily_z)
fft_z[7:(length(fft_z)-6)]<-0
invfftz<-fft(fft_z,inverse=TRUE)/length(fft_z)

lines(day_in_year,daily_z)
lines(day_in_year,invfftz,lwd=2,col="blue")
lines(day_in_year,invfftz+invfft_sd_z,col="purple",lwd=2)
lines(day_in_year,invfftz+150,col="darkgreen",lwd=2)


fft_z2<-fft(daily_z)
fft_z2[21:(length(fft_z2)-20)]<-0
invfftz2<-fft(fft_z2,inverse=TRUE)/length(fft_z2)

#lines(day_in_year,invfftz2,lwd=2,col="purple")

abline(v=59,lwd=2)
abline(v=151,lwd=2)
abline(v=243,lwd=2)
abline(v=334,lwd=2)

text((59+151)/2,5000,"MAM")
text((243+151)/2,5000,"JJA")
text((243+334)/2,5000,"SON")
text(59/2,5000,"DJF")
text((334+366)/2,5000,"DJF")

plot(day_in_year,sd_daily_z,main="Z, standard deviation threshold (daily=black, FFT=blue)\n with constant threshold (green)",
     type="l",xlab="Day in year",ylab="Z (m)")
lines(day_in_year,invfft_sd_z,col="blue",lwd=2)
abline(h=150,lwd=2,col="darkgreen")

abline(v=59,lwd=2)
abline(v=151,lwd=2)
abline(v=243,lwd=2)
abline(v=334,lwd=2)

text((59+151)/2,180,"MAM")
text((243+151)/2,180,"JJA")
text((243+334)/2,180,"SON")
text(59/2,180,"DJF")
text((334+366)/2,180,"DJF")

