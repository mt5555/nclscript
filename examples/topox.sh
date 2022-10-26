#!/bin/tcsh  
if ($#argv < 1) then
  #echo "climoplot.sh casename"
  #exit
endif

set exe = ~/codes/nclscript/contour/vertlevels.py
set wdir = ~/scratch1/preqx/dcmip_tests/dcmip2012_test2.0_steady_state_with_orography/theta-l/movies
set IC = ~/inputdata/atm/cam/inic/homme/cami_mam3_0000-10-01_ne30np4_L72_c160127.nc

#set topo = ~/scratch1/topodata/ne30/USGS-gtopo30_ne30np4pg2_x6t-SGH.nc
#set var=PHIS_d
#set name=x6t

#set topo = ~/scratch1/topodata/ne30/USGS-gtopo30_ne30np4pg2_fx3t-SGH.nc
#set var=PHIS_d
#set name=fx3t

#set topo = ~/scratch1/topodata/ne30/USGS-gtopo30_ne30np4pg2_px4t-SGH.nc
#set var=PHIS_d
#set name=px4t

#set topo = ~/scratch1/topodata/ne30/USGS-gtopo30_ne30pg4-np4pg2_0xdel2-nolim-SGH.nc
#set var=PHIS_d
#set name=x0

set topo = ~/scratch1/topodata/ne30/USGS-gtopo30_ne30np4_16xdel2-nolim.nc
set var=PHIS
set name=x16


conda activate nglmpl2

#$exe -m andes  -i $wdir/nonhydro_dcmip2012_test2_01.nc   geos
#$exe -m andes  -i $topo  -j $IC  $var
$exe -m himalaya  -i $topo  -j $IC  $var
mv ph.png ph-$name.png
mv zh.png zh-$name.png



