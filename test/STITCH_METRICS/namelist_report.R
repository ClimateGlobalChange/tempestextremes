#This is the namelist for the report
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
icfiles<-c("~/BLOBSTATS_FILES/DJF_NA_ic_ERA_MERRA.RData",
           "~/BLOBSTATS_FILES/DJF_NA_ic_ERA_CFSR.RData")#,
           #"~/BLOBSTATS_FILES/DJF_NA_ic_ERA_JRA.RData")


