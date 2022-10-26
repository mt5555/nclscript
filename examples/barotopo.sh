#!/bin/tcsh  
if ($#argv < 1) then
  #echo "climoplot.sh casename"
  #exit
endif

#conda activate nglmpl2
ln -fs  /lcrc/group/acme/data/inputdata/atm/cam/inic/homme/northamericax4v1.g /tmp/temp.nc

cd ~/scratch1/preqx/dcmip_tests/dcmip2016_test1_baroclinic_wave/theta-l

#contour.py -c -.0002,.0002  -m latlon-nc2 -t 7 -i movies/r25-topo_dcmip2016_test11.nc zeta
#\mv -f movies/r25-topo_dcmip2016_test11.zeta.pdf movies/r25-topo_dcmip2016_test11.global.zeta.pdf

#contour.py  -c 730,800,10 -m barotopo2   -t 7  -i movies/r25-topo_dcmip2016_test11.nc ps
#contour.py  -c 940,1020,10  -m latlon-nc2 -t 7 -i movies/r100-topo_dcmip2016_test11.nc \
#-e /lcrc/group/acme/data/inputdata/atm/cam/inic/homme/northamericax4v1.g ps

#/lcrc/group/acme/data/inputdata/atm/cam/inic/homme/northamericax4v1.gcontour.py -c 20,140,1 -m latlon-nc2  -i movies/rrm-topo_dcmip2016_test11.nc ave_dx
#contour.py  -c 940,1020,10  -m barotopo -t 7 -i movies/r25-topo_dcmip2016_test11.nc \
#-e /lcrc/group/acme/data/inputdata/atm/cam/inic/homme/northamericax4v1.g   \
#ps 



#set carg = "-.00025,.00025,.00005"
#contour.py -c $carg -m barotopo -p 750 -t 7 -i movies/r100-topo_dcmip2016_test11.nc zeta
#contour.py -c $carg -m barotopo -p 750 -t 7 -i movies/r25-topo_dcmip2016_test11.nc zeta
#contour.py -c $carg  -m barotopo -p 750 -t 7 -i movies/rrm-topo_dcmip2016_test11.nc zeta

set carg = "230,290,5"
set pval = 700
contour.py  -c $carg -m barotopo -p $pval -t 7 -i movies/r100-topo_dcmip2016_test11.nc T
contour.py  -c $carg -m barotopo -p $pval -t 7 -i movies/r25-topo_dcmip2016_test11.nc T
contour.py  -c $carg -m barotopo -p $pval -t 7 -i movies/rrm-topo_dcmip2016_test11.nc T

#set carg = "230,265,5"
#set pval = 500
#contour.py  -c $carg -m barotopo -p $pval -t 7 -i movies/r100-topo_dcmip2016_test11.nc T
#contour.py  -c $carg -m barotopo -p $pval -t 7 -i movies/r25-topo_dcmip2016_test11.nc T
#contour.py  -c $carg -m barotopo -p $pval -t 7 -i movies/rrm-topo_dcmip2016_test11.nc T
#contour.py  -c $carg -m barotopo -p $pval -t 7 -i movies/nh-rrm-topo_dcmip2016_test11.nc T



#rsync -aP movies/*pdf dosadi:
#rsync -aP movies/*pdf macpro:



