sectors<-c("NA","NC","NP","SA","SI","SP")
seasons<-c("DJF","MAM","JJA","SON")

for (sec in 1:6){
  s<-sectors[sec]
  img_name<-sprintf("~/block_r_data/%s_pv_z_data.RData",s)
  load(img_name)
  ms<-format(time_format,"%b")
  pct90pv<-c()
  pct90z<-c()

  for (v in c("PV","Z")){
    if (v=="PV"){
      xlim=c(0,3e-6)
      ylim=c(0,1.5e6)
    }else{
      xlim=c(0,300)
      ylim=c(0,0.015)
    }
    #Initialize the plot
    png(sprintf("~/Dropbox/plot_anomalies/%s_%s_anom_distr.png",s,v),width=800,height=600)
    plot(NULL,xlim=xlim,ylim=ylim,xlab=sprintf("%s anom",v),ylab="",
         main=sprintf("%s anomalies for %s",v,s))
    for (d in seasons){
      if(d=="DJF"){
        ind_seasons<-which(ms=="Dec" | ms=="Jan" | ms=="Feb")
        lcol<-"blue"
      }else if (d=="MAM"){
        ind_seasons<-which(ms=="Mar" | ms=="Apr" | ms=="May")
        lcol<-"green"
      }else if (d=="JJA"){
        ind_seasons<-which(ms=="Jun" | ms=="Jul" | ms=="Aug") 
        lcol<-"red"
      }else{
        ind_seasons<-which(ms=="Sep"|ms=="Oct"|ms=="Nov")
        lcol<-"purple"
      }
    
      if (v=="PV"){
        var_sub<-pv_anom[,,ind_seasons]
        
        #Only interested in looking at negative (positive) anomalies in NH (SH)
        #interested in magnitude, so make negative anomalies positive
        if (sec<4){
          var_sign<- -var_sub[which(var_sub<0)]
        }else{
          var_sign<-var_sub[which(var_sub>0)]
        }
        abline(v=1.3e-6,lwd=4,lty=2)
      }else{
        var_sub<-z_anom[,,ind_seasons]    
        #Only interested in positive anomalies
        var_sign<-var_sub[which(var_sub>0)]
        abline(v=170,lwd=4,lty=2)
      }
      lines(density(var_sign),col=lcol,lwd=3)
      abline(v=quantile(var_sign,0.9),lty=2,lwd=2,col=lcol)
    }
    dev.off()
  }  
}