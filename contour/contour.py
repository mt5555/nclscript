#!/usr/bin/env python
# Contour over maps of lat/lon data
# python contour_latlon.py
#
#
# contourngl -i inputfile  [-s scripfile]  varname
#
from __future__ import print_function
import os, numpy
import Nio, Ngl
from plotutils import mpl_plot, ngl_plot, myargs
from matplotlib import pyplot



klev=17   # plot level 17 i
proj="latlon"
use_ngl=True

inname,varnames,use_ngl,scrip_file,gll_file = myargs(os.sys.argv)

var1 = varnames[0]
print('contour ',var1,' from ',inname)
infile = Nio.open_file(inname,"r")
dataf  = infile.variables[var1]
print("rank=",dataf.rank," shape=",dataf.shape," dims: ",dataf.dimensions)

title=var1
longname=""
units=""
if hasattr(dataf,"long_name"):
    longname=dataf.long_name
    title=""
if hasattr(dataf,"units"):
    units=dataf.units
        
    
if "ncol_d" in dataf.dimensions:
    lat  = infile.variables["lat_d"][:]
    lon  = infile.variables["lon_d"][:]
else:
    lat  = infile.variables["lat"][:]
    lon  = infile.variables["lon"][:]



# does array have a time dimension?
ntimes=1
times=numpy.array([0])
if "time" in dataf.dimensions:
    ntimes = infile.dimensions['time']
    times = infile.variables["time"][:]
    

if use_ngl:
    wks_type = "png"
    wks = Ngl.open_wks(wks_type,"testplot")
    cmap="notyetcodedforngl"
else:
    cmap='nipy_spectral'


clev=numpy.array([50])   # 50 levels, no range specified
clev=numpy.array([-.0005,.0005,.0005/20])


for t in range(ntimes):
    timedim =  "time" in dataf.dimensions
    levdim = "lev" in dataf.dimensions

    if (timedim and levdim):
        data2d=dataf[t,klev,...]
        print("time=",times[t]," k=",klev)
    elif (timedim):
        data2d=dataf[t,...]
        print("time=",times[t])
    elif (levdim):
        data2d=dataf[klev,...]
        print("k=",klev)
    else:
        data2d=dataf.values

    if use_ngl:
        ngl_plot(wks,wks_type,data2d,lon,lat,title,longname,units,
                 proj,clev,cmap,scrip_file)
    else:
        mpl_plot(data2d,lon,lat,title,longname,units,
                 proj,clev,cmap,gll_file)
        pyplot.savefig("testplot.png",dpi=300)
        #pyplot.savefig("testplot.pdf",dpi=300)
        #pyplot.show()

if use_ngl:
    del wks
    Ngl.end()
