# This file generates correlation coefficients, Q-Q plots, and histogram
# plots for the dataframe generated in sector.R
# Pre
#Need to run this first if data not run:
# source("~/tempestextremes/test/sector.R")

#Break methods: "FD", "even_number", or "even_space"

chart_stats<-function(varname,data_vec,labs,seasons_vec,sec,hemi,
                      colvec=c("blueviolet","blue","green","aquamarine2","cyan4"), 
                      calc_corr=FALSE,qq=FALSE,hist_plot=FALSE,cdf_plot=FALSE,box_plot=FALSE,
                      printb=FALSE,ks_test=FALSE,dval=FALSE,pval=FALSE, xsq=FALSE,ad=FALSE,
                      breaks_method="even_number",break_interval=2, num_breaks=40,save_img=FALSE){
  #Check to see if the proper directory exists
  img_dir=paste("~/Desktop/Research reports/",varname,"/",sep="")
  if (!dir.exists(img_dir)){
    dir.create(img_dir)
  }
  
  #Setting up the various parameters
  nYearCompare=20
  nPerCompare=4
  
  ########################
  
  # Entering the loop 
  short=c("ERA","Cl","2C","S2","CS")
  for (season in seasons_vec){
    for (h in hemi){
      for(s in sec){
        if (calc_corr==TRUE){
          print(paste("Spearman correlation between datasets for",varname,season,h,s))
        }
        dat=df_tot[df_tot$Sector==s & df_tot$Hemi==h & df_tot$Season==season,
                   c("centlon_c","Dataset","nYears",varname)]

        if (varname=="centlon"){
          if (s=="PAC"){
            varname_plot="centlon"
          }else if (s=="ATL"){
            varname_plot="centlon_c"
          }
        }else if (varname != "centlon"){
          varname_plot=varname
        }
        climo=dat[dat$Dataset=='climo',varname_plot]
        ERA=dat[dat$Dataset=='ERA',varname_plot]
        SST=dat[dat$Dataset=='SSTplus2',varname_plot]
        CO2=dat[dat$Dataset=='2xCO2',varname_plot]
        SST_CO2=dat[dat$Dataset=='SSTplus2_2xCO2',varname_plot]
        if(ad==TRUE){
          dats=list(climo,ERA,SST,CO2,SST_CO2)
          for (x in dats){
            for (y in dats){
              adt=ad.test(x,y)
            }
          }
        }
        
        if (printb==TRUE){
          #b=adjbox(climo,ERA,CO2,SST,SST_CO2,names=data_vec,plot=FALSE)
          b=boxplot(climo,ERA,CO2,SST,SST_CO2,names=data_vec,plot=FALSE)
          b_stats=data.frame(b$stats, row.names=c("lw","lh","med","uh","uw"))
          colnames(b_stats)=data_vec               
          return(b_stats)
          #Each row of b$stats: lower whisker, lower hinge, median, upper hinge, extreme            
        }
        if (box_plot==TRUE){
          if (save_img==TRUE){
            pname=paste(img_dir,season,"_",h,"_",s,"_",varname,"_boxplot.png",sep="")
            png(pname,width=1200,height=900)
          }
          #print(data_vec)
          mt=paste(varname,"boxplots for",season,h,s)
          b=adjbox(climo,ERA,CO2,SST,SST_CO2,names=short,las=2, 
                   main=mt,cex.labs=1.8,cex.main=1.3,col=colvec)
          if (save_img==TRUE){
            dev.off()
          }

        }
        else{
          #Get the extent of values for the breaks vector
          lower=floor(min(dat[,varname_plot]))
          upper=ceiling(max(dat[,varname_plot]))
  
          #Generate the vector of breaks
          if (breaks_method=="FD"){
            nCounts<-c()
            interq<-c()
            for (d in 1:length(data_vec)){
              sub=dat[dat$Dataset==data_vec[d],]
              if (data_vec[d]=='ERA'){
                t=4
              }else{
                t=8
              }
              n=sub$nYears[1]
              scale=(nYearCompare*nPerCompare)/(t*n)
              nPoints=nrow(sub)
              nCounts<-c(nCounts,scale*nPoints)
              interq<-c(interq,IQR(sub[,varname_plot]))
            }
            nAvg=mean(nCounts)
            iqAvg=mean(interq)
  
            break_interval=round(2*iqAvg/(nAvg^(1/3)),1)
          }
          else if (breaks_method=="even_number"){
            break_interval=round((upper-lower)/num_breaks, 1)
          }
          else if (breaks_method=="even_space"){
            break_interval=break_interval
          }
  
          breaks=seq(lower,upper+break_interval,break_interval)
  
          mat=data.frame(matrix(nrow=length(breaks)-1,ncol=length(data_vec)))
          colnames(mat)<-data_vec
          
          #Get the histogram max for ylim
          maxh=0
          for (d in 1:length(data_vec)){
            sub=dat[dat$Dataset==data_vec[d],] 
            h1=hist(sub[,varname_plot],breaks=breaks,plot=FALSE,include.lowest=TRUE)
  
            maxh1=max(h1$density)
            maxh=max(maxh,maxh1)
          }
          
          
          for (d in 1:length(data_vec)){
            sub=dat[dat$Dataset==data_vec[d],]
            if (data_vec[d]=="ERA"){
              t=4
            }else{
              t=8
            }
  
            h1=hist(sub[,varname_plot],breaks=breaks,plot=FALSE,include.lowest=TRUE)
            mat[,d]=h1$density
            if(xsq==TRUE){
              for (x in c(1,3:5)){
                xtest=chisq.test(mat[,2],mat[,x])
              }
            }
            if (hist_plot==TRUE){
              if (save_img==TRUE){
                pname=paste(img_dir,season,"_",h,"_",s,"_",varname,
                            "_",breaks_method,"_", data_vec[d],"_hist.png",sep="")
                png(pname,width=1200,height=900)
              }
              xlim=range(breaks)
              mt=paste(data_vec[d],"blocking",varname,"density for",season,h,s)
              break_label=paste("Histogram break:",break_interval)
              ylim=c(0,maxh)
              plot(h1,main=mt,sub=break_label,xlab=varname,col=colvec[d],
                   xlim=xlim,ylim=ylim,freq=FALSE,cex.main=1.3)
              if (save_img==TRUE){
                dev.off()
              }
              if (data_vec[d]=='ERA'){
                #print(breaks)
                if (breaks_method=="FD"){
                  print(nCounts)
                  print(interq)
                  print(c(nAvg,iqAvg))
                }
              }
            }
          }
          xvec=c(min(na.omit(mat)),max(na.omit(mat)))
          if (ks_test==TRUE){
            ksmat_d=data.frame(matrix(nrow=length(data_vec),ncol=length(data_vec)),row.names=labs)
            colnames(ksmat_d)<-labs
            ksmat_p=data.frame(matrix(nrow=length(data_vec),ncol=length(data_vec)),row.names=labs)
            colnames(ksmat_p)<-labs
            for(i in 1:length(data_vec)){
              for (j in 1:length(data_vec)){
                k=ks.test(mat[,i],mat[,j])
                ksmat_d[i,j]=k$statistic
                ksmat_p[i,j]=round(k$p.value,4)
              }
            }
            if (dval==TRUE){
              return(ksmat_d)
            }
            if (pval==TRUE){
              return(ksmat_p)
            }
  
          }
          if (cdf_plot==TRUE){
            if (save_img==TRUE){
              pname=paste(img_dir,season,"_",h,"_",s,"_",varname,"_",breaks_method,"_cdf_compare.png",sep="")
              png(pname,width=1200,height=900)
            }
  
            par(mar=c(4,4,1,1)+.1)
            mt=paste("Blocking",varname,"CDFs for",season,h,s)
            plot(ecdf(mat[,1]),main=mt,col=colvec[1],xlab='density',cex.lab=1.6)
            for (j in 2:length(data_vec)){
              lines(ecdf(mat[,j]),main=mt,col=colvec[j],cex.lab=1.6)
            }
            legend(x=(xvec[2]-0.5*(xvec[2]-xvec[1])),y=0.2,legend=labs,col=colvec,lwd=3)
            if (save_img==TRUE){
              dev.off()
            }
          }
          if (qq==TRUE){
            if (save_img==TRUE){
              pname=paste(img_dir,season,"_",h,"_",s,"_",varname, "_", breaks_method,"_climo_compare.png",sep="")
              png(pname,width=1200,height=300)
            }
  
            par(mfrow=c(1,4))
            par(mar=c(4,4,1,1)+.1)
            for (i in 2:2){
              for (j in 1:length(data_vec)){
                if (i!=j){
                  qqplot(mat[,j],mat[,i],xlim=xvec,ylim=xvec,xlab=labs[j],ylab=labs[i],cex.lab=1.6)
                  abline(a=0,b=1)
                }
              }
            }
            if (save_img==TRUE){
              dev.off()
            }
          }
          if (calc_corr==TRUE){
            print(cor(mat,method="pearson"))
          }
        }
      }
    }
  }
}

