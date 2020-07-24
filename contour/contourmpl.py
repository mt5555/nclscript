#!/usr/bin/env python
# Contour over maps of lat/lon data
# python contour_latlon.py
#
# needs my nglutils module:
# setenv PYTHONPATH ${HOME}/codes/nclscript/modules
#
# contourngl -i inputfile  [-g gllfile]  varname
#
from __future__ import print_function
import numpy,os,sys,getopt
import Nio
from matplotlib import pyplot
from plotutils import myargs, mpl_plot



inname,varnames,scripfile,gllfile = myargs(os.sys.argv)

var1 = varnames[0]
print('contour ',var1,' from ',inname)
infile = Nio.open_file(inname,"r")

data2d  = infile.variables[var1]
lat  = infile.variables["lat"][:]
lon  = infile.variables["lon"][:]

# set contours: vmin,vmax,nlevels

# colormap:
cmap='nipy_spectral'

# projection:
proj="latlon"
#proj="andes"
clev=numpy.array([50]);  # 50 levels, no range specified
mpl_plot(data2d,lon,lat,"title",proj,clev,cmap,gllfile)

#pyplot.savefig("testplot.png",dpi=300)
#pyplot.savefig("testplot.pdf",dpi=300)
pyplot.show()

