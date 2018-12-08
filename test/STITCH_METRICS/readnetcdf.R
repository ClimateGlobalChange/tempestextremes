require(ncdf4)
require(ncdf4.helpers)
require(RNetCDF)
require(abind)
require(akima)
require(PCICt)

#Returns in the -180->180 range
lon_convert<-function(lon){
  distFrom180=lon-180.
  return(ifelse(
    distFrom180<0,
    lon,
    -(180-distFrom180)
  ))
}
#Returns in the 0->360 range
lon_convert2<-function(lon){
  return(ifelse(lon<0,360+lon,lon))
}

#Finds the index of closest value in a vector
nearest_ind<-function(val,vec){
  dist_vec<-abs(val-vec)
  i<-which(dist_vec==min(dist_vec))
  return(i[1])
}


linint<-function(from_arr,from_lat,from_lon,to_lat,to_lon){
  to_lon<-round(to_lon,3)
  to_lat<-round(to_lat,3)
  from_lon<-round(from_lon,3)
  from_lat<-round(from_lat,3)
  dx<-to_lon[2]-to_lon[1]
  dy<-to_lat[2]-to_lat[1]
  # print(sprintf("%f %f",dx,dy))
  new_grid<-bilinear.grid(x=from_lon,y=from_lat,z=from_arr,
                          xlim=range(to_lon),ylim=range(to_lat),
                          dx=dx,dy=dy)
  return(new_grid$z)
}


read_netcdf<-function(flist,vlist,olist=vlist,timename="time",levname="lev",latname="lat",lonname="lon",
                      minlat="",maxlat="",minlon="",maxlon="",minlev="",maxlev="",ncout="",rdataout="",
                      transformto180=FALSE,transformto360=FALSE,regridto1degree=FALSE){
  
  minlat<-ifelse(minlat=="",NA,as.numeric(minlat))
  
  maxlat<-ifelse(maxlat=="",NA,as.numeric(maxlat))
  minlon<-ifelse(minlon=="",NA,as.numeric(minlon))
  maxlon<-ifelse(maxlon=="",NA,as.numeric(maxlon))
  minlev<-ifelse(minlev=="",NA,as.numeric(minlev))
  maxlev<-ifelse(maxlev=="",NA,as.numeric(maxlev))
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
    print(sprintf("Reading file %s",f))
    ncfile<-nc_open(f)
    num_dims<-ncfile$ndims
    #print(sprintf("file has %d dimensions",num_dims))
    #Get the names of each dimension
    #Note that the dimension id indexing starts at 0
    dimnames<-nc.get.dim.names(ncfile)
    #Create a variable for each of the available dimensions
    in_dims<-c(timename,levname,latname,lonname)
    lat_units<-ncatt_get(ncfile,latname,"units")$value
    lon_units<-ncatt_get(ncfile,lonname,"units")$value
    
    assigned_names<-c("timevar","lev","lat","lon")
    axis_names<-c("time_axis","lev_axis","lat_axis","lon_axis")
    axis_type<-c()
    axis_units<-c()
    saved_names<-c()
    for (i in 1:4){
      if (sum(dimnames==in_dims[i])>0){
        #Create a variable with a name that matches the NetCDF variable
        #vartype<-var.inq.nc(ncfile,in_dims[i])$type
        vartype<-"NC_DOUBLE"
        assign(assigned_names[i],ncvar_get(ncfile,in_dims[i]))
        saved_names<-c(saved_names,axis_names[i])
        axis_type<-c(axis_type,vartype)
        axis_units<-c(axis_units,ncfile$dim[[in_dims[i]]]$units)
      }
    }
    time_units<-ncfile$dim[[timename]]$units
    calstring<-ncatt_get(ncfile,timename,"calendar")$value
    print(sprintf("calendar is %s",calstring))
    #print(time_units)
    #print(timevar[1])
    tstring<-as.character(nc.get.time.series(ncfile))
    #tstring<-utcal.nc(time_units,timevar,type="s")
    #print(tstring[1])
    #time_axis<-c(time_axis,timevar)
    time_format<-c(time_format,tstring)
    #What is the direction of the latitude axis?
    POSTOP<-FALSE
    if (lat[1]-lat[2]>0){
      POSTOP<-TRUE
    }
    if (transformto180==TRUE){
      lon<-lon_convert(lon)
    }
    if (transformto360==TRUE){
      lon<-lon_convert2(lon)
    }
    #if (subset==TRUE){
    #			if (is.null(minlat) | is.null(maxlat) | is.null(minlon) | is.null(maxlon)){
    #stop("Need to specify lat/lon boundaries for subset option.")
    #			}
    if (is.na(minlat) & is.na(maxlat)){
      bi<-1
      ti<-length(lat)
    }else{
      bi<-nearest_ind(minlat,lat)
      ti<-nearest_ind(maxlat,lat)
    }
    if (length(bi)<1 | length(ti)<1){
      stop("Can't find specified min and max latitude coordinates.")
    }
    lat_inds<-seq(bi,ti)
    #print(lat_inds)
    lat_axis<-lat[lat_inds]
    #print(lat_axis)
    
    if (is.na(minlon) & is.na(maxlon)){
      li<-1
      ri<-length(lon)
    }else{
      li<-nearest_ind(minlon,lon)
      ri<-nearest_ind(maxlon,lon)
      #print(sprintf("Subsetting bounds to %f, %f",lon[li],lon[ri]))
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
      if (!is.na(minlev) & !is.na(maxlev)){
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
    
    v_units<-c()
    for (i in 1:length(vlist)){
      #Read in these variables
      #Check the dimensionality
      vcheck<-ncfile$var[[vlist[i]]]
      v_units<-c(v_units,vcheck$units)
      #dimids 
      dims_present<-dimnames[(vcheck$dimids)+1]
      ndims<-vcheck$ndims
      v_in<-NULL
      if (ndims==3){
        #lon,lat,time
        v_in<-ncvar_get(ncfile,vlist[i])[lon_inds,lat_inds,]	
      }else if (ndims==4){
        #lon,lat,lev,time
        v_in<-ncvar_get(ncfile,vlist[i])[lon_inds,lat_inds,lev_inds,]
      }
      #concatenate along time axis (last axis)
      vbind<-abind(get(olist[i]),v_in)
      assign(olist[i],vbind)			
    }	
    nc_close(ncfile)
  }
  if (is.unsorted(time_format)){
    t_sortorder<-sort(time_format,index.return=T)$ix
    tcopy<-time_axis

    for (v in olist){
      ndims<-length(dim(get(v)))
      if (ndims==3){
        v_sort<-get(v)[,,t_sortorder]
      }else if (ndims==4){
        v_sort<-get(v)[,,,t_sortorder]
      }
      assign(v,v_sort)
    }

    time_format<-tstringcopy[t_sortorder]
  }
  time_axis<-utinvcal.nc("hours since 1800-01-01 00:00",time_format)

  if (regridto1degree==TRUE){
    vlist_orig<-sprintf("%s_orig",olist)
    #Save the original dimension information-- time will be the same
    lat_axis_orig<-lat_axis
    lon_axis_orig<-lon_axis
    lat_axis<-seq(as.integer(minlat),as.integer(maxlat))
    lon_axis<-seq(as.integer(minlon),as.integer(maxlon))
    
    for (i in 1:length(olist)){
      #Copy the variable to a new name
      temp_var<-get(olist[i])
      print(sprintf("regridding %s",olist[i]))
      assign(vlist_orig[i],temp_var)
      ndims<-length(dim(temp_var))
      #Regrid the data
      new_var<-array(0,dim=c(length(lon_axis),length(lat_axis),length(time_format)))
      for (t in 1:length(time_format)){
        if (ndims==3){
          calc_var<-temp_var[,,t]
          new_var[,,t]<-linint(calc_var,lat_axis_orig,lon_axis_orig,lat_axis,lon_axis)
        }else{
          #print("Will only regrid for 3 dimension variables.")
        }
      }
      new_var[which(new_var>0)]<-1
      #new_var[which(new_var<0.5)]<-0
      assign(olist[i],new_var)
    }
    olist_fin<-c(olist,vlist_orig,"lat_axis_orig","lon_axis_orig")
  }else{
    olist_fin<-olist
  }
  save_variables<-c("time_format",saved_names,olist_fin)
  if (ncout!= ""){
    dimlist<-list()
    ncvarlist<-list()
    for (i in 1:length(saved_names)){
      dimlen<-length(get(saved_names[i]))
      tempdim<-ncdim_def(saved_names[i],axis_units[i],get(saved_names[i]),calendar = "standard")
      dimlist[[i]]<-tempdim
    }

    
    for (i in 1:length(olist)){
      dimorder<-ncfile$var[[vlist[i]]]$dimids
      dim_copy<-dimlist[1:length(dimorder)]
      for (n in 1:length(dimorder)){
        dim_copy[[dimorder[n]+1]]<-dimlist[[n]]
      }
      tempvar<-ncvar_def(name=olist[i],units=v_units[i],dim=dim_copy)
      ncvarlist[[i]]<-tempvar
      
    }
    ncfile_out<-nc_create(ncout,ncvarlist)
    #ncatt_put(ncfile_out,"time_axis","calendar","standard")
    for (i in 1:length(olist)){
      ncvar_put(ncfile_out,olist[i],get(olist[i]))
    }
    nc_close(ncfile_out)
    print(sprintf("Wrote %s to file",ncout))
  }
  if (rdataout!=""){
    save(list=save_variables,file=rdataout)
    print(sprintf("Wrote %s to file",rdataout))
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