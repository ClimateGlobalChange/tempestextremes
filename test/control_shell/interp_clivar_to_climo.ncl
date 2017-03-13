load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin
	in=addfile(fin,"r")
	lon = in->lon(:)
	lat = in->lat(:)
	
	newlon=fspan(0.,359.,360)
	newlon@units=lon@units
	newlat=fspan(-90.,90.,181)
	newlat@units=lat@units
	
	;;Read in the variables
	newipv=linint2(lon,in->IPV&lat,in->IPV,True,newlon,newlat,0)
	newipv!0="lat"
	newipv!1="lon"
	newipv&lat=newlat
	newipv&lon=newlon
	
	newT=linint2(lon,in->AVGT&lat,in->AVGT,True,newlon,newlat,0)
	newT!0="lat"
	newT!1="lon"
	newT&lat=newlat
	newT&lon=newlon
	
	newU=linint2(lon,in->AVGU&lat,in->AVGU,True,newlon,newlat,0)
	newU!0="lat"
	newU!1="lon"
	newU&lat=newlat
	newU&lon=newlon
	
	newV=linint2(lon,in->AVGV&lat,in->AVGV,True,newlon,newlat,0)
	newV!0="lat"
	newV!1="lon"
	newV&lat=newlat
	newV&lon=newlon
	
	system("/bin/rm " + fout)
	out=addfile(fout,"c")
	
	out->IPV = newipv
	out->AVGT = newT
	out->AVGU = newU
	out->AVGV = newV
	
end