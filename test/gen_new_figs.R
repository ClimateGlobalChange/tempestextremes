library(RNetCDF)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(reshape2)
library(RColorBrewer)
library(xtable)
library(scales)
source("~/tempestextremes/test/sector_funcs.R")


brks_z<-seq(4500,6100,50)
hgt.cols<-colorRampPalette(c("purple","blue","cyan4","green","yellow", "orange","red","darkred"))(length(brks_z))

#JJA STUFF----------
load('~/block_r_data/ERA_NA_pv_z_ghg_block_data.RData')

time_hr1<-time_format
time_format1<-substr(time_hr1,1,10)
load('~/block_r_data/ERA_NA_z_inst_data.RData')
time_hr2<-time_format

time_hrs<-rep(seq(0,18,6),length(time_hr1)/4)
#Need to reverse the latitude axis for each of the data
lats_seq<-lats_seq[length(lats_seq):1]
z_inst<-z_inst[,length(lats_seq):1,]
z_anom<-z_anom[,length(lats_seq):1,]
pv_anom<-pv_anom[,length(lats_seq):1,]
ghg<-ghg[,length(lats_seq):1,]


lons_plot=lon_convert(lons_seq)
mapr="world"
m<-map_data(mapr)


ratio.display <- 4/3
ratio.values <- (max(lons_seq)-min(lons_seq))/(max(lats_seq)-min(lats_seq))
subfigs<-c("a","b","c","d")

ctr<-1
print(range(lons_plot))
dates_plot<-c("1984-06-08_12","1984-06-10_12","1984-06-12_12","1984-06-14_12")
for (x in dates_plot){
  t1<-which(time_hr1==x)
  t2<-which(time_hr2==x)
  longdata_sub<-melt(z_inst[,,t2],varnames=c("x","y"))
  longdata_sub$lon<-lons_plot[longdata_sub[,1]]
  longdata_sub$lat<-lats_seq[longdata_sub[,2]]
  longdata_sub$cont<-melt(ghg[,,t1])$value
  longdata_sub$cont2<-melt(pv_anom[,,t1])$value
  longdata_sub$cont3<-melt(z_anom[,,t1])$value
  
  g<-ggplot()+ 
    coord_fixed(ratio.values / ratio.display) +
    coord_cartesian(xlim=c(min(lons_plot),max(lons_plot)),
                    ylim=c(min(lats_seq),max(lats_seq)),expand=FALSE) +
    geom_map(data= m, map = m, aes(map_id=region)) +
    geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=value,color=..level..),breaks=brks_z,size=1) +
    scale_color_gradientn(limits=c(min(brks_z),max(brks_z)),colors=hgt.cols,name="Z500 (m)") +
    geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=cont),breaks=c(0,1),color="purple",size=3.25) +
    geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=cont2),breaks=c(0,1),color="cornflowerblue",size=3.25) +
    geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=cont3),breaks=c(0,1),color="chartreuse4",size=3.25) +
    geom_rect(aes(xmin = -105, xmax=-50,ymin=30,ymax=45),
              fill = "transparent", color = "blue", size = 1.5) +
    ggtitle(sprintf("(%s) %s %sZ",subfigs[ctr],time_format1[t1],time_hrs[t1])) +
    labs(x="Longitude",y="Latitude")+
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title=element_text(size=22),
          legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(2,"line"),
          legend.key.width=unit(3.5,"line"),
          plot.title=element_text(size=24)) 
  #png(sprintf("~/PAPERS_WRITING/PAPER1/lolat_JJA/NA_%sZ_noT.png",time_hr1[t]),width=800,height=450)
  png(sprintf("~/block_r_data/figs/fig14%s.png",subfigs[ctr]),width=600,height=450)
  print(g)
  dev.off()
  ctr<-ctr+1
}
#Also need to do wind and temperature data--------
#load the netcdfs for zonal data
#
 zonbrk<-seq(-20,20,2)
 zon_cols<-colorRampPalette(c("green","turquoise4","darkblue",
                              "white","darkred","purple","orange"))(length(zonbrk))
#
#
 zonalnc<-open.nc("~/block_r_data/ERA_1984_06_vars_U_devs.nc")
 lat<-read.nc(zonalnc)$lat
 lon<-read.nc(zonalnc)$lon
 time_str<-as.Date(read.nc(zonalnc)$time/24,origin="1800-01-01")
 hrs_vec<-rep(seq(0,18,6),(length(time_str)/4))
 time_hrstr<-sprintf("%s_%02d",time_str,hrs_vec)
 t1<-which(time_hrstr=="1984-06-08_12")
 t2<-which(time_hrstr=="1984-06-17_00")
 zon3d<-read.nc(zonalnc)$ADU[,,t1:t2]
 close.nc(zonalnc)

 zon<-apply(zon3d,c(1,2),mean)
 dat<-melt(zon)
 dat$lon<-lon[dat[,1]]
 dat$lon2<-lon_convert(dat$lon)
 dat$lat<-lat[dat[,2]]
#
 meridnc<-open.nc("~/block_r_data/ERA_1984_06_vars_V_devs.nc")
 mer3d<-read.nc(meridnc)$ADV[,,t1:t2]
 close.nc(meridnc)
 merid<-apply(mer3d,c(1,2),mean)
 dat$ADV<-melt(merid)$value
 dat$brks<-breaks_cuts(dat,zonbrk)
#
 mapr="world"
 m<-map_data(mapr)
#
 gz<-ggplot()+
   #  coord_fixed(ratio.values / ratio.display) +
   coord_cartesian(xlim=c(-110,-45),
                   ylim=c(25,50),expand=FALSE) +
   geom_map(data= m, map = m, aes(map_id=region)) +
   geom_raster(data=dat,aes(x=lon2,y=lat,fill=value),interpolate = T,alpha=0.5) +
   geom_contour(data=dat,aes(x=lon2,y=lat,z=value,color=..level..),breaks=zonbrk,size=1) +
   scale_fill_gradientn(breaks=zonbrk,limits=c(-20,20),
                        colors=zon_cols,name="U anom (m/s)",
                        labels=c("-20",rep("",9),"0",rep("",9),20)) +
   scale_color_gradientn(breaks=zonbrk,limits=c(-20,20),colors=zon_cols,guide=F) +
   geom_rect(aes(xmin = -105, xmax=-50,ymin=30,ymax=45),
             fill = "transparent", color = "blue", size = 1.5) +
   ggtitle(" (c) Average zonal wind anomaly\n1984-06-06 to 1984-06-17") +
   labs(x="Longitude",y="Latitude")+
   theme(axis.text=element_text(size=28),
         axis.title=element_text(size=30),
         legend.text=element_text(size=28),
         legend.title=element_text(size=30),
         legend.direction = "horizontal",
         legend.position = "bottom",
         legend.key.size = unit(2,"line"),
         legend.key.width=unit(3.5,"line"),
         plot.title=element_text(size=32))

# #png("~/PAPERS_WRITING/PAPER1/lolat_JJA/zonal_anom_JJA.png",width=600,height=450)
 png("~/PAPERS_WRITING/PAPER1_REV/fig11c.png",width=600,height=500)
 print(gz)
 dev.off()
#
 merbrk<-seq(-16,16,2)
 mer_cols<-colorRampPalette(c("orange","darkgreen",
                              "white","purple","violetred4"))(length(merbrk))
#
gm<-ggplot()+
  #  coord_fixed(ratio.values / ratio.display) +
  coord_cartesian(xlim=c(-110,-45),
                  ylim=c(25,50),expand=FALSE) +
  geom_map(data= m, map = m, aes(map_id=region)) +
  geom_raster(data=dat,aes(x=lon2,y=lat,fill=ADV),interpolate = T,alpha=0.5) +
  geom_contour(data=dat,aes(x=lon2,y=lat,z=ADV,color=..level..),breaks=merbrk,size=1) +
  scale_fill_gradientn(breaks=merbrk,limits=c(-16,16),
                       colors=mer_cols,name="V anom (m/s)",
                       labels=c("-16",rep("",7),"0",rep("",7),16)) +
  scale_color_gradientn(breaks=merbrk,limits=c(-20,20),colors=mer_cols,guide=F) +
  geom_rect(aes(xmin = -105, xmax=-50,ymin=30,ymax=45),
            fill = "transparent", color = "blue", size = 1.5) +
  ggtitle(" (b) Average meridional wind anomaly\n1984-06-06 to 1984-06-17") +
  labs(x="Longitude",y="Latitude")+
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30),
        legend.text=element_text(size=28),
        legend.title=element_text(size=30),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.size = unit(2,"line"),
        legend.key.width=unit(3.5,"line"),
        plot.title=element_text(size=32))

#png("~/PAPERS_WRITING/PAPER1/lolat_JJA/zonal_anom_JJA.png",width=600,height=450)
png("~/PAPERS_WRITING/PAPER1_REV/fig11b.png",width=600,height=500)
print(gm)
dev.off()




#DJF STUFF-------------


#Also need to do wind and temperature data--------
#load the netcdfs for zonal data
#
zonbrk<-seq(-20,20,2)
zon_cols<-colorRampPalette(c("green","turquoise4","darkblue",
                             "white","darkred","purple","orange"))(length(zonbrk))


zonalnc<-open.nc("~/block_r_data/ERA_2006_01_vars_U_devs.nc")
lat<-read.nc(zonalnc)$lat
lon<-read.nc(zonalnc)$lon
time_str<-as.Date(read.nc(zonalnc)$time/24,origin="1800-01-01")
hrs_vec<-rep(seq(0,18,6),(length(time_str)/4))
time_hrstr<-sprintf("%s_%02d",time_str,hrs_vec)
t1<-which(time_hrstr=="2006-01-05_00")
t2<-which(time_hrstr=="2006-01-11_00")
zon3d<-read.nc(zonalnc)$ADU[,,t1:t2]
close.nc(zonalnc)

zon<-apply(zon3d,c(1,2),mean)
dat<-melt(zon)
dat$lon<-lon[dat[,1]]
dat$lat<-lat[dat[,2]]

meridnc<-open.nc("~/block_r_data/ERA_2006_01_vars_V_devs.nc")
mer3d<-read.nc(meridnc)$ADV[,,t1:t2]
close.nc(meridnc)
merid<-apply(mer3d,c(1,2),mean)
dat$ADV<-melt(merid)$value
#dat$brks<-breaks_cuts(dat,zonbrk)

mapr="world2"
m<-map_data(mapr)

gz<-ggplot()+ 
  #  coord_fixed(ratio.values / ratio.display) +
  coord_cartesian(xlim=c(150,230),
                  ylim=c(25,50),expand=FALSE) +
  geom_map(data= m, map = m, aes(map_id=region)) +
  geom_raster(data=dat,aes(x=lon,y=lat,fill=value),interpolate = T,alpha=0.5) +
  geom_contour(data=dat,aes(x=lon,y=lat,z=value,color=..level..),breaks=zonbrk,size=1) +
  scale_fill_gradientn(breaks=zonbrk,limits=c(-20,20),
                       colors=zon_cols,name="U anom (m/s)",
                       labels=c("-20",rep("",9),"0",rep("",9),"20"),oob=squish) +
  scale_color_gradientn(breaks=zonbrk,limits=c(-20,20),colors=zon_cols,guide=F) +
  geom_rect(aes(xmin = 160, xmax=220,ymin=26,ymax=47),
            fill = "transparent", color = "blue", size = 1.5) +
  ggtitle("(c) Average zonal wind anomaly\n2006-01-05 to 2006-01-11") +
  labs(x="Longitude",y="Latitude")+
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30),
        legend.text=element_text(size=28),
        legend.title=element_text(size=20),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.size = unit(2,"line"),
        legend.key.width=unit(3.5,"line"),
        plot.title=element_text(size=32)) 

#png("~/PAPERS_WRITING/PAPER1/lolat_JJA/zonal_anom_JJA.png",width=600,height=450)
png("~/PAPERS_WRITING/PAPER1_REV/fig13c.png",width=600,height=450)
print(gz)
dev.off()
#
merbrk<-seq(-16,16,2)
mer_cols<-colorRampPalette(c("orange","darkgreen",
                             "white","purple","violetred4"))(length(merbrk))
#
gm<-ggplot()+ 
  #  coord_fixed(ratio.values / ratio.display) +
  coord_cartesian(xlim=c(150,230),
                  ylim=c(25,50),expand=FALSE) +
  geom_map(data= m, map = m, aes(map_id=region)) +
  geom_raster(data=dat,aes(x=lon,y=lat,fill=ADV),interpolate = T,alpha=0.5) +
  geom_contour(data=dat,aes(x=lon,y=lat,z=ADV,color=..level..),breaks=merbrk,size=1) +
  scale_fill_gradientn(breaks=merbrk,limits=c(-16,16),
                       colors=mer_cols,name="V anom (m/s)",
                       labels=c("-16",rep("",7),"0",rep("",7),"16")) +
  scale_color_gradientn(breaks=merbrk,limits=c(-16,16),colors=mer_cols,guide=F) +
  geom_rect(aes(xmin = 160, xmax=220,ymin=26,ymax=47),
            fill = "transparent", color = "blue", size = 1.5) +
  ggtitle("(b) Average merid wind anomaly\n2006-01-05 to 2006-01-11") +
  labs(x="Longitude",y="Latitude")+
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=30),
        legend.text=element_text(size=28),
        legend.title=element_text(size=30),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.size = unit(2,"line"),
        legend.key.width=unit(3.5,"line"),
        plot.title=element_text(size=32)) 

# #png("~/PAPERS_WRITING/PAPER1/lolat_JJA/zonal_anom_JJA.png",width=600,height=450)
png("~/PAPERS_WRITING/PAPER1_REV/fig13b.png",width=600,height=450)
print(gm)
dev.off()






#Attempt at wind plots

windbrk<-seq(0,50,5)
wind_cols<-colorRampPalette(c("white","blue","green","yellow","orange","red"))(length(windbrk))

#JJA--------
u<-open.nc("~/block_r_data/ERA_1984_06_vars.nc")
lev<-read.nc(u)$lev
p<-which(lev==500)
time_str<-as.Date(read.nc(u)$time/24,origin="1800-01-01")
hrs_vec<-rep(seq(0,18,6),(length(time_str)/4))
time_hrstr<-sprintf("%s_%02d",time_str,hrs_vec)
t1<-which(time_hrstr=="1984-06-08_12")
t2<-which(time_hrstr=="1984-06-17_00")
umat<-read.nc(u)$U[,,p,t1:t2]
vmat<-read.nc(u)$V[,,p,t1:t2]
lat<-read.nc(u)$lat
lon<-read.nc(u)$lon
close.nc(u)

uw<-apply(umat,c(1,2),mean)
vw<-apply(vmat,c(1,2),mean)


windtab<-melt(uw)
windtab$lat<-lat[windtab[,2]]
windtab$lon<-lon[windtab[,1]]
windtab$lon2<-lon_convert(windtab$lon)
windtab$v<-melt(vw)$value

windtab$magnitude<-sqrt(windtab$value^2 + windtab$v^2)
windtab$theta<-atan(windtab$v/windtab$value)
windtab$theta2<-ifelse(windtab$value<0,windtab$theta+pi,windtab$theta)
windtab$ulength<-cos(windtab$theta2)
windtab$vlength<-sin(windtab$theta2)


gw<-ggplot(windtab,aes(x=lon2,y=lat, xend=lon2+1.5*ulength,yend=lat+1.5*vlength)) + 
  coord_cartesian(xlim=c(-110,-45),
                  ylim=c(25,50),expand=FALSE) +
  geom_raster(aes(fill=magnitude),interpolate = T)+
  scale_fill_gradientn(breaks=windbrk,limits=c(0,50),
                       colors=wind_cols,name="Total wind 500 mb (m/s)") +
  geom_segment(data=windtab[windtab$lat%%5==0 & windtab$lon%%5==0,],
               arrow=arrow(length = unit(0.05, "cm"),type="closed") )+
  geom_rect(aes(xmin = -105, xmax=-50,ymin=30,ymax=45),
            fill = "transparent", color = "blue", size = 1.5) +
  ggtitle("(a) Total wind 500 mb,  1984-06-06 to 1984-06-17")+
  labs(x="Longitude",y="Latitude")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title=element_text(size=22),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.size = unit(2,"line"),
        legend.key.width=unit(3.5,"line"),
        plot.title=element_text(size=24)) 

#png("~/PAPERS_WRITING/PAPER1/lolat_JJA/wind_JJA.png",width=600,height=450)
png("~/PAPERS_WRITING/PAPER1_REV/fig18a.png",width=600,height=450)
print(gw)
dev.off()


#DJF----
u<-open.nc("~/block_r_data/ERA_2006_01_vars.nc")
time_str<-as.Date(read.nc(u)$time/24,origin="1800-01-01")
hrs_vec<-rep(seq(0,18,6),(length(time_str)/4))
time_hrstr<-sprintf("%s_%02d",time_str,hrs_vec)
t1<-which(time_hrstr=="2006-01-05_00")
t2<-which(time_hrstr=="2006-01-11_00")
umat<-read.nc(u)$U[,,p,t1:t2]
vmat<-read.nc(u)$V[,,p,t1:t2]
lat<-read.nc(u)$lat
lon<-read.nc(u)$lon
close.nc(u)

uw<-apply(umat,c(1,2),mean)
vw<-apply(vmat,c(1,2),mean)

windtab<-melt(uw)

windtab$lat<-lat[windtab[,2]]
windtab$lon<-lon[windtab[,1]]
windtab$v<-melt(vw)$value

windtab$magnitude<-sqrt(windtab$value^2 + windtab$v^2)
windtab$theta<-atan(windtab$v/windtab$value)
windtab$theta2<-ifelse(windtab$value<0,windtab$theta+pi,windtab$theta)
windtab$ulength<-cos(windtab$theta2)
windtab$vlength<-sin(windtab$theta2)

gw<-ggplot(windtab,aes(x=lon,y=lat, xend=lon+1.5*ulength,yend=lat+1.5*vlength)) + 
  coord_cartesian(xlim=c(150,230),
                  ylim=c(25,50),expand=FALSE) +
  geom_raster(aes(fill=magnitude),interpolate = T)+
  scale_fill_gradientn(breaks=windbrk,limits=c(0,50),
                       colors=wind_cols,name="Total wind 500 mb (m/s)") +
  geom_segment(data=windtab[windtab$lat%%5==0 & windtab$lon%%5==0,],
               arrow=arrow(length = unit(0.05, "cm"),type="closed") )+
  geom_rect(aes(xmin = 160, xmax=220,ymin=26,ymax=47),
            fill = "transparent", color = "blue", size = 1.5) +
  ggtitle("(b) Total wind 500 mb, 2006-01-05 to 2006-01-11")+
  labs(x="Longitude",y="Latitude")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text=element_text(size=20),
        legend.title=element_text(size=22),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.key.size = unit(2,"line"),
        legend.key.width=unit(3.5,"line"),
        plot.title=element_text(size=24)) 

#png("~/PAPERS_WRITING/PAPER1/lolat_DJF/wind_DJF.png",width=600,height=450)
png("~/PAPERS_WRITING/PAPER1_REV/fig18b.png",width=600,height=450)
print(gw)
dev.off()

#WIND SHEAR-----------


load('~/block_r_data/old_data/NA_pv_z_ghg_block_data.RData')

time_format1<-time_format
time_hr1<-sprintf("%s_%02d",time_format1,time_hrs)
lats_seq<-lats_seq[length(lats_seq):1]
pv_anom<-pv_anom[,length(lats_seq):1,]

lons_plot=lons_seq_c
mapr="world"
m<-map_data(mapr)

dates_plot<-c("1989-09-28_06","1989-09-29_06","1989-09-30_06","1989-10-01_06")

#make a matrix to hold the data
umat<-array(NA,c(360,181,4))
vmat<-array(NA,c(360,181,4))
#Read the first three dates

u<-open.nc("~/block_r_data/ERA_1989_09_vars.nc")
time_str<-as.Date(read.nc(u)$time/24,origin="1800-01-01")
hrs_vec<-rep(seq(0,18,6),(length(time_str)/4))
time_hrstr<-sprintf("%s_%02d",time_str,hrs_vec)
t1<-which(time_hrstr==dates_plot[1])
t2<-which(time_hrstr==dates_plot[2])
t3<-which(time_hrstr==dates_plot[3])
lat<-read.nc(u)$lat
lon<-read.nc(u)$lon
umat[,,1:3]<-read.nc(u)$U[,,p,c(t1,t2,t3)]
vmat[,,1:3]<-read.nc(u)$V[,,p,c(t1,t2,t3)]
close.nc(u)

u<-open.nc("~/block_r_data/ERA_1989_10_vars.nc")
time_str<-as.Date(read.nc(u)$time/24,origin="1800-01-01")
hrs_vec<-rep(seq(0,18,6),(length(time_str)/4))
time_hrstr<-sprintf("%s_%02d",time_str,hrs_vec)
t1<-which(time_hrstr==dates_plot[4])
umat[,,4]<-read.nc(u)$U[,,p,c(t1)]
vmat[,,4]<-read.nc(u)$V[,,p,c(t1)]
close.nc(u)

nth<-c("1st","2nd","3rd","4th")
fn<-c("sep28","sep29","sep30","oct1")
for (x in 1:4){
  
  t<-which(time_hr1==dates_plot[x])
  longdata_sub<-melt(pv_anom[,,t],varnames=c("x","y"))
  longdata_sub$lon<-lons_plot[longdata_sub[,1]]
  longdata_sub$lat<-lats_seq[longdata_sub[,2]]
  

  uw<-umat[,,x]
  vw<-vmat[,,x]
  windtab<-melt(uw)

  
  windtab$lat<-lat[windtab[,2]]
  windtab$lon<-lon[windtab[,1]]
  windtab$lon2<-ifelse(windtab$lon>180,windtab$lon-360,windtab$lon)
  windtab$v<-melt(vw)$value
  
  windtab$magnitude<-sqrt(windtab$value^2 + windtab$v^2)
  windtab$theta<-atan(windtab$v/windtab$value)
  windtab$theta2<-ifelse(windtab$value<0,windtab$theta+pi,windtab$theta)
  windtab$ulength<-cos(windtab$theta2)
  windtab$vlength<-sin(windtab$theta2)
  
  
  gw<-ggplot() + 
    coord_cartesian(xlim=c(-100,40),
                    ylim=c(25,75),expand=FALSE) +
    geom_raster(data=windtab,aes(x=lon2,y=lat,fill=magnitude),interpolate = T)+
    scale_fill_gradientn(breaks=windbrk,limits=c(0,50),
                         colors=wind_cols,name="Total wind 500 mb (m/s)",oob=squish) +
    geom_segment(data=windtab[windtab$lat%%5==0 & windtab$lon%%5==0,],
                 aes(x=lon2,y=lat, xend=lon2+1.5*ulength,yend=lat+1.5*vlength),
                 arrow=arrow(length = unit(0.05, "cm")) )+
    geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=value),
                 breaks=c(0,1),color="black",size=3.25) +
    #geom_rect(aes(xmin = -105, xmax=-50,ymin=30,ymax=45),
    #          fill = "transparent", color = "blue", size = 1.5) +
    ggtitle(sprintf("(%s) Total wind 500 mb,  %s",subfigs[x],substr(dates_plot[x],1,10)))+
    labs(x="Longitude",y="Latitude")+
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=22),
          legend.text=element_text(size=20),
          legend.title=element_text(size=22),
          legend.direction = "horizontal",
          legend.position = "bottom",
          legend.key.size = unit(2,"line"),
          legend.key.width=unit(3.5,"line"),
          plot.title=element_text(size=24)) 
  
  #png(sprintf("~/PAPERS_WRITING/PAPER1/pv_diff/%s.png",fn[x]),width=600,height=450)
   png(sprintf("~/PAPERS_WRITING/PAPER1_REV/fig11%s.png",subfigs[x]),width=600,height=450)
  print(gw)
   dev.off()
}
