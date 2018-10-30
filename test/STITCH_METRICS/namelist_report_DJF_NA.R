#This is the namelist for the report

metadata_datasets="This is the results summary for DJF NA (December 1980-February 2010). 
JRA is 1.25x1.25 degree, ERA is 1x1 degree, MERRA is 0.5x0.625 degree, and CFSR is 0.5x0.5 degree"

#This will be the output name
output_name<-"~/BLOBSTATS_FILES/DJF_NA_allreanalysis.html"
#Load the datasets
Varnames<-c("JRA","ERA","MERRA","CFSR")
#Number of years in the dataset
nyears<-2009-1980+1
#Per-timestep data
mergefiles<-sprintf("~/BLOBSTATS_FILES/%s/%s_DJF_NA_merged_table.RData",Varnames,Varnames)
#Summary data
summfiles<-sprintf("~/BLOBSTATS_FILES/%s/%s_DJF_NA_summ_table.RData",Varnames,Varnames)
#Blob files
blobfiles<-sprintf("~/BLOBSTATS_FILES/%s/%s_DJF_NA_blobdata.RData",Varnames,Varnames)
#Name of the blob variable
blobname<-rep("Z_BLOB",4)

#Intercomparison data
#Ideally, do it in the order (1,2) (1,3) (1,4) (2,3) (2,4) (3,4) so that the table is nicely organized!
icfiles<-c()
for (x in 1:length(Varnames)){
  for (y in 1:length(Varnames)){
    if (x<y){
      fn<-sprintf("~/BLOBSTATS_FILES/DJF_NA_ic_%s_%s.RData",Varnames[x],Varnames[y])
      icfiles<-c(icfiles,fn)
    }
  }
}
# 
# icfiles<-c("~/BLOBSTATS_FILES/DJF_NA_ic_JRA_ERA.RData",
#           "~/BLOBSTATS_FILES/DJF_NA_ic_ERA_MERRA.RData",
#            "~/BLOBSTATS_FILES/DJF_NA_ic_ERA_CFSR.RData")#,
#            #"~/BLOBSTATS_FILES/DJF_NA_ic_ERA_JRA.RData")


