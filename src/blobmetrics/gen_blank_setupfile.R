#This script will generate a blank namelist which can then be filled in

fname<-"blank_setupfile.R"
sink(fname,split=T)

cat("#The requirements are as follows:
    #1) StitchBlobs files from each dataset
    #2) Corresponding BlobStats files from each dataset
    #3) Optional but recommended: BlobStats output from 
    # the separate DetectBlobs function 
    # (This will help deal with cases where two blobs 
    # merge into a single blob at a later time)  
    
    #REPORT HEADER: ANY TEXT THAT YOU WANT TO GO AT THE TOP OF THE REPORT 
    #NOTE: must be in quotes. Careful of any special characters 
    metadata_datasets=\"\" 
    #This is the number of years in your dataset! (Integer) 
    nyears<-0  

    ########### 
    #FILE INFO 
    #Will always output RData files 
    #But there is an option to also output text files 
    # or CSV files with the table data 
    #Set following variables to TRUE if you wish to have one or both of these outputs 
    output_txt<-FALSE 
    output_csv<-FALSE  

    #Names of the datasets (vector of strings) 
    #Example: Varnames<-c(\"ERA-Interim\",\"CFSR\",\"MERRA-2\") 
    #For all of the vectors of strings, make sure that the lengths are identical to 
    # the length of this vector! 
    Varnames<-c()  

    #The spatial resolutions of the datasets (Vector of strings.  
    #Optional-- delete this variable if you don't want to use it) 
    #Example: resolutions<-c(\"1x1\",\"0.5x0.625\",\"0.5x0.5\") 
    resolutions<-c()  

    #Working directory (String) 
    #This is where all of the R function files are stored  
    work_dir<-\"\" 
    #Output directory(String) 
    # Main directory where all output files will go
    output_dir<-\"\" 
    #Name of output subdirectory (String) 
    # This will be created in the output directory specified above 
    output_subdir<-\"\" 
    #Name of the file prefix (String).  
    #Output file will be [prefix]_[Varname]_[suffix] (for example, [prefix]_ERA_stitchtable.RData)  
    output_prefix<-\"\" 
    #Output name for the master namelist that will be generated using this file template (String) 
    fname_namelist<-\"\"  


    #Input lists of files 
    #Must use full pathname for each of these lists! 
    #Example: stitch_lists<-c(\"~/input_dir/ERA/ERA_list\",\"~/input_dir/CFSR/CFSR_list\",\"~/input_dir/MERRA/MERRA_list\") 
    #List of BlobStats files (from StitchBlobs output) to read into readfiles.R (vector of strings) 
    stitch_lists<-c() 
    #Use DetectBlobs inputs? (TRUE/FALSE) 
    use_detectblob<-TRUE  
    #List of BlobStats files (from DetectBlobs output) to read into readfiles.R (vector of strings) 
    #If not using BlobStats files from DetectBlobs, then delete this variable 
    detect_lists<-c() 
    #List of StitchBlobs files to read into readnetcdf.R (vector of strings) 
    stitchblob_lists<-c() 
    #Name of the StitchBlobs variable in the NetCDF file (vector of strings) 
    #Note: if the variable name is identical for all file, you can do 
    # varvec<-rep(\"VARNAME\",length(Varnames)) 
    varvec<-c() 
    #Name of the time, lat, lon axes (String) 
    timename<-\"time\" 
    latname<-\"lat\" 
    lonname<-\"lon\" 
    #Transform the lon axis? (TRUE/FALSE) 
    #Note: if a dataset already has this longitude extent, it will do nothing 
    #from 0/360 to -180/180 
    transformto180<-FALSE 
    #from -180/180 to 0/360 
    transformto360<-TRUE 
    #Subset lat and lon if desired (note: minlon and maxlon should correspond to 
    # the appropriate longitude extent, i.e. either in the -180/180 range or 0/360 range) 
    #If you don't wish to subset, delete these four variables 
    minlat<- -90 
    maxlat<- 90 
    minlon<- 0 
    maxlon<- 359 
    #Regrid to 1 degree? (TRUE/FALSE) 
    regridto1degree<-TRUE
    ")
sink()
