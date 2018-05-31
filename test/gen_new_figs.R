library(RNetCDF)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(reshape2)
library(RColorBrewer)
library(xtable)
source("~/tempestextremes/test/sector_funcs.R")


brks_z<-seq(4500,6100,50)
hgt.cols<-colorRampPalette(c("purple","blue","cyan4","green","yellow", "orange","red","darkred"))(length(brks_z))

load('~/block_r_data/ERA_NP_pv_z_ghg_block_data.RData')

time_hr1<-time_format
time_format1<-substr(time_hr1,1,10)
load('~/block_r_data/ERA_NP_z_inst_data.RData')
time_hr2<-time_format

time_hrs<-rep(seq(0,18,6),length(time_hr1)/4)
#Need to reverse the latitude axis for each of the data
lats_seq<-lats_seq[length(lats_seq):1]
z_inst<-z_inst[,length(lats_seq):1,]
z_anom<-z_anom[,length(lats_seq):1,]
pv_anom<-pv_anom[,length(lats_seq):1,]
ghg<-ghg[,length(lats_seq):1,]

lons_plot=lons_seq
mapr="world2"
m<-map_data(mapr)


ratio.display <- 4/3
ratio.values <- (max(lons_seq)-min(lons_seq))/(max(lats_seq)-min(lats_seq))
subfigs<-c("a","b","c","d")

ctr<-1
print(range(lons_plot))
#Figure lowlat blocking, with new boundaries
dates_plot<-c("1991-07-29_12","1991-08-01_12","1991-08-04_12","1991-08-07_12")
for (x in dates_plot){
  t1<-which(time_hr1==x)
  t2<-which(time_hr2==x)
  longdata_sub<-melt(z_inst[,,t2],varnames=c("x","y"))
  longdata_sub$lon<-lons_plot[longdata_sub[,1]]
  longdata_sub$lat<-lats_seq[longdata_sub[,2]]
  longdata_sub$cont1<-melt(z_anom[,,t1])$value
  longdata_sub$cont2<-melt(pv_anom[,,t1])$value
  longdata_sub$cont3<-melt(ghg[,,t1])$value

  g<-ggplot()+
    coord_fixed(ratio.values / ratio.display) +
    coord_cartesian(xlim=range(lons_plot),
                    ylim=range(lats_seq),expand=FALSE) +
    geom_map(data= m, map = m, aes(map_id=region)) +
    geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=value,color=..level..),breaks=brks_z,size=1) +
    scale_color_gradientn(limits=c(min(brks_z),max(brks_z)),colors=hgt.cols,name="Z500 (m)") +
    geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=cont1),breaks=c(0,1),color="cornflowerblue",size=3.25) +
    geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=cont3),breaks=c(0,1),color="purple",size=3.25) +
    geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=cont2),breaks=c(0,1),color="chartreuse4",size=3.25) +
    geom_rect(aes(xmin = 153, xmax=220,ymin=32,ymax=48),
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
          plot.title=element_text(size=24,hjust=0))
  #png(sprintf("~/PAPERS_WRITING/PAPER1/z_zg/SP_%sZ_noT.png",time_hr1[t]),width=800,height=450)
 # png(sprintf("~/PAPERS_WRITING/CLIM_DYN_PAPER1/fig12%s.png",subfigs[ctr]),width=600,height=450)
  png(sprintf("~/block_r_data/figs/fig14%s.png",subfigs[ctr]),width=600,height=450)
  print(g)
  dev.off()
  ctr<-ctr+1
}

#Also need to do wind and temperature data--------
#load the netcdfs for zonal data
zonbrk<-seq(-10,10,2)
zon_cols<-colorRampPalette(c("darkblue","white","darkred"))(length(zonbrk))


zonalnc<-open.nc("/group/paullricgrp4/ERA_MCP/ERA_1991/wind_devs.nc")
lat<-read.nc(zonalnc)$lat
lon<-read.nc(zonalnc)$lon
time<-read.nc(zonalnc)$time
time_string<-as.Date(time/24, origin="1800-01-01")
time_hours<-rep(seq(0,18,6),(length(time)/4))
time_strhr<-sprintf("%s_%02d",time_string,time_hours)
t1<-which(time_strhr=="1991-07-29_12")
t2<-which(time_strhr=="1991-08-09_00")
zondev<-read.nc(zonalnc)$ADU[,,t1:t2]
meriddev<-read.nc(zonalnc)$ADV[,,t1:t2]
close.nc(zonalnc)

zon<-apply(zondev,c(1,2),mean)
mer<-apply(meriddev,c(1,2),mean)
dat<-melt(zon,value.name = "ADU")
dat$lon<-lon[dat[,1]]
dat$lon2<-ifelse(dat$lon>180,dat$lon-360,dat$lon)
dat$lat<-lat[dat[,2]]
dat$ADV<-melt(mer)$value


mapr="world"
m<-map_data(mapr)

gz<-ggplot()+ 
  #  coord_fixed(ratio.values / ratio.display) +
  coord_cartesian(xlim=c(-75,-15),
                  ylim=c(25,55),expand=FALSE) +
  geom_map(data= m, map = m, aes(map_id=region)) +
  geom_raster(data=dat,aes(x=lon2,y=lat,fill=ADU),interpolate = T,alpha=0.5) +
  geom_contour(data=dat,aes(x=lon2,y=lat,z=ADU,color=..level..),breaks=zonbrk,size=1) +
  scale_fill_gradientn(breaks=zonbrk,limits=c(-10,10),
                       colors=zon_cols,name="U anom (m/s)",
                       labels=c("-10",rep("",4),"0",rep("",4),10)) +
  scale_color_gradientn(breaks=zonbrk,limits=c(-10,10),colors=zon_cols,guide=F) +
  geom_rect(aes(xmin = -70, xmax=-22,ymin=32,ymax=50),
            fill = "transparent", color = "blue", size = 1.5) +
  ggtitle(" (c) Average zonal wind anomaly\n2007-07-29 to 2007-08-04") +
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

#png("~/PAPERS_WRITING/PAPER1/lolat_JJA/zonal_anom_JJA.png",width=600,height=450)
png("~/block_r_data/figs/fig15c.png",width=600,height=450)
print(gz)
dev.off()


#Figure RRR
# ctr<-1
# dates_plot<-c("1989-09-28_06","1989-09-29_06","1989-09-30_06","1989-10-01_06")
# for (x in dates_plot){
#   t1<-which(time_hr1==x)
#   t2<-which(time_hr2==x)
#   longdata_sub<-melt(z_inst[,,t2],varnames=c("x","y"))
#   longdata_sub$lon<-lons_plot[longdata_sub[,1]]
#   longdata_sub$lat<-lats_seq[longdata_sub[,2]]
#   longdata_sub$cont1<-melt(z_anom[,,t1])$value
#   longdata_sub$cont2<-melt(pv_anom[,,t1])$value
#   longdata_sub$cont3<-melt(ghg[,,t1])$value
# 
#   g<-ggplot()+
#     coord_fixed(ratio.values / ratio.display) +
#     coord_cartesian(xlim=c(min(lons_plot),max(lons_plot)),
#                     ylim=c(min(lats_seq),max(lats_seq)),expand=FALSE) +
#     geom_map(data= m, map = m, aes(map_id=region)) +
#     geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=value,color=..level..),breaks=brks_z,size=1) +
#     scale_color_gradientn(limits=c(min(brks_z),max(brks_z)),colors=hgt.cols,name="Z500 (m)") +
#     geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=cont1),breaks=c(0,1),color="cornflowerblue",size=3.25) +
#     geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=cont3),breaks=c(0,1),color="purple",size=3.25) +
#     geom_contour(data=longdata_sub,aes(x=lon,y=lat,z=cont2),breaks=c(0,1),color="chartreuse4",size=3.25) +
#     ggtitle(sprintf("(%s) %s %sZ",subfigs[ctr],time_format1[t1],time_hrs[t1])) +
#     labs(x="Longitude",y="Latitude")+
#     theme(axis.text=element_text(size=20),
#           axis.title=element_text(size=22),
#           legend.text=element_text(size=20),
#           legend.title=element_text(size=22),
#           legend.direction = "horizontal",
#           legend.position = "bottom",
#           legend.key.size = unit(2,"line"),
#           legend.key.width=unit(3.5,"line"),
#           plot.title=element_text(size=24,hjust=0))
#   #png(sprintf("~/PAPERS_WRITING/PAPER1/z_zg/SP_%sZ_noT.png",time_hr1[t]),width=800,height=450)
#  # png(sprintf("~/PAPERS_WRITING/CLIM_DYN_PAPER1/fig12%s.png",subfigs[ctr]),width=600,height=450)
#   png(sprintf("~/block_r_data/figs/fig10%s.png",subfigs[ctr]),width=600,height=450)
#   print(g)
#   dev.off()
#   ctr<-ctr+1
# }
# 
# 
# 
