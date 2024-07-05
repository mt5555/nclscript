#!/bin/tcsh
#
# script for plotting the soluation of the no-perturbation JW instability
# test, i.e. from LJTN James 2009 paper
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.3894/JAMES.2010.2.15
#
# in LJTN paper, used 997,1003,.5 contour intervals at day 9
#
contour.py -i steady-tensor.nc -c 997,1003,.5 -t 3 -y mpl -f la -o 0 ps
contour.py -i steady-nu1e15.nc -c 997,1003,.5 -t 3 -y mpl -f la -o 0 ps
contour.py -i steady-nu5e14.nc -c 997,1003,.5 -t 3 -y mpl -f la -o 0  ps
rename .ps.pdf -- -t9.ps.pdf steady-*.ps.pdf 


contour.py -i steady-tensor.nc -c 999.5,1000.5,.1 -y mpl -f la -o 0 ps
contour.py -i steady-nu1e15.nc -c 999.5,1000.5,.1 -y mpl -f la -o 0 ps
contour.py -i steady-nu5e14.nc -c 999.5,1000.5,.1 -y mpl -f la -o 0 ps

# NGL version:
#contour.py -i steady-tensor.nc -c 999.5,1000.5,.1  -f la -o 0 ps
#contour.py -i steady-nu1e15.nc -c 999.5,1000.5,.1  -f la -o 0 ps
#contour.py -i steady-nu5e14.nc -c 999.5,1000.5,.1  -f la -o 0 ps

