#This is the namelist for the report

metadata_datasets="This is the results summary for JJA SP (June 1980-August 2009). 
JRA is 1.25x1.25 degree, ERA is 1x1 degree, MERRA is 0.5x0.625 degree, and CFSR is 0.5x0.5 degree"

output_name<-"~/BLOBSTATS_FILES/JJA_SP_allreanalysis.html"
#Load the datasets
Varnames<-c("JRA","ERA","MERRA","CFSR")
#Number of years in the dataset
nyears<-2009-1980+1
#Per-timestep data
mergefiles<-sprintf("~/BLOBSTATS_FILES/%s/%s_JJA_SP_merged_table.RData",Varnames,Varnames)
#Summary data
summfiles<-sprintf("~/BLOBSTATS_FILES/%s/%s_JJA_SP_summ_table.RData",Varnames,Varnames)
#Blob files
blobfiles<-sprintf("~/BLOBSTATS_FILES/%s/%s_JJA_SP_blobdata.RData",Varnames,Varnames)
#Name of the blob variable
blobname<-rep("Z_BLOB",4)

#Intercomparison data
#Ideally, do it in the order (1,2) (1,3) (1,4) (2,3) (2,4) (3,4) so that the table is nicely organized!
icfiles<-c()
for (x in 1:length(Varnames)){
  for (y in 1:length(Varnames)){
    if (x<y){
      fn<-sprintf("~/BLOBSTATS_FILES/JJA_SP_ic_%s_%s.RData",Varnames[x],Varnames[y])
      icfiles<-c(icfiles,fn)
    }
  }
}

# icfiles<-c("~/BLOBSTATS_FILES/JJA_SP_ic_ERA_MERRA.RData",
#            "~/BLOBSTATS_FILES/JJA_SP_ic_ERA_CFSR.RData")#,
#            #"~/BLOBSTATS_FILES/JJA_SP_ic_ERA_JRA.RData")


