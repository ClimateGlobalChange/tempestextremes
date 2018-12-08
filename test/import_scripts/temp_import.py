from ecmwfapi import ECMWFDataServer
import sys
 
server = ECMWFDataServer()
 
server.retrieve({
 "dataset"   : "interim",
    "date"      : sys.argv[1], #dates
    "levelist"  : "250/500/850", #includes information on these pressure levels
    "stream"    : "oper",
    "format"    : "netcdf",
    "levtype"   : "pl",
    "param"     : "129.128/130.128/131.128/132.128", #PV, T, U, V, rel vorticity
    "step"      : "0",
    "time"      : "00:00:00/06:00:00/12:00:00/18:00:00", 
    "type"      : "an",
    "area"      : "90/0/-90/359",
    "grid"      : "1.0/1.0",
    "target"    : sys.argv[2],
   })
sys.exit(0)
