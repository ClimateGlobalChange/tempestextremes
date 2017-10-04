#!/bin/bash
let a=0
for f in IVT_MERRA_201702/*.nc
do
  ../bin/SpineARs --in $f --out output/$f --var IVT --minlaplacian 50000 --laplaciansize 10 --minval 250 --minabslat 15 --zonalmeanwt 0 --meridmeanwt 0 --minarea 40 --addtimedim $a --addtimedimunits "minutes since 2017-02-01 00:00:00"
  let a+=180
done

ncrcat output/IVT_MERRA_201702/* MERRA2.ar_tag.Tempest_v1.3hourly.20170201-20170228.nc
