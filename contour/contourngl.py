#!/usr/bin/env python
# Contour over maps of lat/lon data
# python contour_latlon.py
#
# needs my nglutils module:
# setenv PYTHONPATH ${HOME}/codes/nclscript/modules
#
# contourngl -i inputfile  [-s scripfile]  varname
#
from __future__ import print_function
import os, numpy
import Nio, Ngl
from plotutils import ngl_plot, myargs

inname,varnames,scrip_file,gll_file = myargs(os.sys.argv)

var1 = varnames[0]
print('contour ',var1,' from ',inname)
infile = Nio.open_file(inname,"r")



data2d  = infile.variables[var1]
lat  = infile.variables["lat"]
lon  = infile.variables["lon"]

proj="latlon"
wks_type = "png"
wks = Ngl.open_wks(wks_type,"testplot")

cmap="notyetcodedforngl"
clev=numpy.array([50])   # 50 levels, no range specified
map=ngl_plot(wks,wks_type,data2d,lon,lat,"Vortex test field",proj,clev,cmap,scrip_file)


del map
del wks
Ngl.end()
