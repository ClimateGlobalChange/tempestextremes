require(RNetCDF)
require(abind)

options(warn=-1)
#This program reads in NetCDFs based on the specified variable list
#The time axis is sorted if out of order
#The output is either 1) RData file with variable names that correspond to the original
#2) A NetCDF file

read_netcdf<-function(flist,vlist,olist=vlist,timename="time",levname="lev",latname="lat",lonname="lon",
	minlat="",maxlat="",minlon="",maxlon="",minlev="",maxlev="",ncout="",rdataout=""){

  minlat<-ifelse(minlat=="",NULL,as.numeric(minlat))
  maxlat<-ifelse(maxlat=="",NULL,as.numeric(maxlat))
  minlon<-ifelse(minlon=="",NULL,as.numeric(minlon))
  maxlon<-ifelse(maxlon=="",NULL,as.numeric(maxlon))
  minlev<-ifelse(minlev=="",NULL,as.numeric(minlev))
  maxlev<-ifelse(maxlev=="",NULL,as.numeric(maxlev))
	#Create a time axis for appending subsequent values
	#time_axis<-c()
	time_format<-c()
	#Create a variable for each of the specified ones in the list
	for (v in olist){
		assign(v,NULL)
	}
	
	#Open each netCDF variable
	#Note that default axis order is lon,lat,lev,time!
	for (f in flist){
		ncfile<-open.nc(f)
		num_dims<-file.inq.nc(ncfile)$ndims

		#Get the names of each dimension
		#Note that the dimension id indexing starts at 0
		dim_info<-sapply(seq(0,(num_dims-1)),dim.inq.nc,ncfile=ncfile)
		dimnames<-unlist(dim_info["name",])
		#Create a variable for each of the available dimensions
		in_dims<-c(timename,levname,latname,lonname)
		lat_units<-att.get.nc(ncfile,latname,"units")
		lon_units<-att.get.nc(ncfile,lonname,"units")
		
		assigned_names<-c("timevar","lev","lat","lon")
		axis_names<-c("time_axis","lev_axis","lat_axis","lon_axis")
		axis_type<-c()
		saved_names<-c()
		for (i in 1:4){
			if (sum(dimnames==in_dims[i])>0){
				#Create a variable with a name that matches the NetCDF variable
				vartype<-var.inq.nc(ncfile,in_dims[i])$type
				assign(assigned_names[i],var.get.nc(ncfile,in_dims[i]))
				saved_names<-c(saved_names,axis_names[i])
				axis_type<-c(axis_type,vartype)
			}
		}
		time_units<-att.get.nc(ncfile,timename,"units")
		tstring<-utcal.nc(time_units,timevar,type="s")
		#time_axis<-c(time_axis,timevar)
		time_format<-c(time_format,tstring)
		#What is the direction of the latitude axis?
		POSTOP<-FALSE
		if (lat[1]-lat[2]>0){
			POSTOP<-TRUE
		}
		
		#if (subset==TRUE){
#			if (is.null(minlat) | is.null(maxlat) | is.null(minlon) | is.null(maxlon)){
				#stop("Need to specify lat/lon boundaries for subset option.")
#			}
		if (is.null(minlat) & is.null(maxlat)){
		  bi<-1
		  ti<-length(lat)
		}else{
		  bi<-which(lat==minlat)
		  ti<-which(lat==maxlat)
		}
		if (length(bi)<1 | length(ti)<1){
		  stop("Can't find specified min and max latitude coordinates.")
		}
		lat_inds<-seq(bi,ti)
		#print(lat_inds)
		lat_axis<-lat[lat_inds]
		#print(lat_axis)
		
		if (is.null(minlon) & is.null(maxlon)){
		  li<-1
		  ri<-length(lon)
		}else{
		  li<-which(lon==minlon)
		  ri<-which(lon==maxlon)
		}
		if (length(li)<1 | length(ri)<1){
		  stop("Can't find specified min and max longitude coordinates.")
		}
		#Do we need to deal with the periodic boundary condition?
		if (li>ri){
		  left_part<-seq(li,length(lon))
		  right_part<-seq(1,ri)
		  lon_inds<-c(left_part,right_part)
		}else{
		  lon_inds<-seq(li,ri)
		}
		lon_axis<-lon[lon_inds]

		#Deal with lev axis, if it exists
		#Check if lev_axis is in the list of saved names
		levpos<-which(saved_names=="lev_axis")
		if (length(levpos)>0){
		  if (!is.null(minlev) & !is.null(maxlev)){
		    #print(sprintf("subsetting level axis from %s to %s",minlev,maxlev))
		    l1<-which(lev==minlev)
		    l2<-which(lev==maxlev)
		    lev_inds<-seq(l1,l2)
		    lev_axis<-lev[lev_inds]
		  }else{
		    lev_inds<-seq(1,length(lev))
		    lev_axis<-lev
		  }
		}

		for (i in 1:length(vlist)){
			#Read in these variables
			#Check the dimensionality
			vcheck<-var.inq.nc(ncfile,vlist[i])
			#dimids 
			dims_present<-dimnames[(vcheck$dimids)+1]
			ndims<-length(dims_present)
			v_in<-NULL
			if (ndims==3){
				#lon,lat,time
				v_in<-var.get.nc(ncfile,vlist[i])[lon_inds,lat_inds,]	
			}else if (ndims==4){
				#lon,lat,lev,time
				v_in<-var.get.nc(ncfile,vlist[i])[lon_inds,lat_inds,lev_inds,]
			}
			#concatenate along time axis (last axis)
			vbind<-abind(get(olist[i]),v_in)
			assign(olist[i],vbind)			
		}	
		close.nc(ncfile)
	}
	#Check to make sure the time steps are in the correct order!!
	#print(is.unsorted(tvec))
	if (is.unsorted(time_format)){
		#print("sorting time axis")
		t_sortorder<-sort(time_format,index.return=T)$ix
		#print("sorted time axis")
		#tcopy<-time_axis
		tstringcopy<-time_format
		#print("copied times")
		for (v in olist){
			ndims<-length(dim(get(v)))
			if (ndims==3){
				v_sort<-get(v)[,,t_sortorder]
			}else if (ndims==4){
				v_sort<-get(v)[,,,t_sortorder]
			}
			assign(v,v_sort)
			#print(dim(get(v)))
		}
		
		#time_axis<-tcopy[t_sortorder]
		time_format<-tstringcopy[t_sortorder]
	}
	#print(time_format)
	time_axis<-utinvcal.nc("hours since 1800-01-01 00:00",time_format)
	#print(time_axis)
	save_variables<-c("time_format",saved_names,olist)
	#print(save_variables)
	
	if (ncout!= ""){
		ncfile_out<-create.nc(ncout)
		#Write the dimensions
		for (i in 1:length(saved_names)){
			dimlen<-length(get(saved_names[i]))
			#print(sprintf("length of variable %s is %d and axis type is %s",saved_names[i],dimlen,axis_type[i]))
			dim.def.nc(ncfile_out,saved_names[i],dimlen)
			var.def.nc(ncfile_out,saved_names[i],axis_type[i],saved_names[i])
			var.put.nc(ncfile_out,saved_names[i],get(saved_names[i]))
		}
		#Define the time units
		att.put.nc(ncfile_out,"time_axis","units","NC_CHAR","hours since 1800-01-01 00:00")
		att.put.nc(ncfile_out,"lat_axis","units","NC_CHAR",lat_units)
		att.put.nc(ncfile_out,"lon_axis","units","NC_CHAR",lon_units)
		#Write the variables to file
		for (i in 1:length(olist)){
			ndims<-length(dim(get(olist[i])))
			if (ndims==3){
				var.def.nc(ncfile_out,olist[i],"NC_DOUBLE",c("lon_axis","lat_axis","time_axis"))
			}else if (ndims==4){
				var.def.nc(ncfile_out,olist[i],"NC_DOUBLE",c("lon_axis","lat_axis","lev_axis","time_axis"))
			}else{
				stop("Currently only writing 3 and 4 dimension variables. Contact developer.")
			}
			#print(sprintf("defined variable %s",vlist[i]))
			var.put.nc(ncfile_out,olist[i],get(olist[i]))

		}
		close.nc(ncfile_out)
	}
	if (rdataout!=""){
		save(list=save_variables,file=rdataout)
	}

	if (rdataout=="" & ncout==""){
		#Create a list to return all of the saved objects
		outlist<-list()
		outlist[["save_variables"]]<-save_variables
		for (v in save_variables){
			outlist[[v]]<-get(v)
		}
		return(outlist)
	}

}