#!/usr/bin/env python
#
# vertlevels.py  -i datafile(hyam,ps) -j PHIS_file   [options]  varname
#
# plot vertical levels
#
#
from __future__ import print_function
import os, numpy
import Nio
from plotutils import myargs, extract_level, interp_to_latlon
from matplotlib import pyplot



inname,inname2,varnames,proj,timeindex,klev,plev,clev,nlatlon_interp,use_ngl,scrip_file,gll_file,se_file \
    = myargs(os.sys.argv)

var2_read="PHIS"
units=""
longname=""
print('vcoord file=',inname)
print('PHIS file=',inname2)
infile = Nio.open_file(inname,"r")
outname=inname2.split(".nc")[0] + ".topo"

nlat=1
nlon=1024
lat_i = numpy.array([-20.0])
lon_i = numpy.linspace(0, 360, nlon,endpoint=False)
#nlat=1024
#nlon=1
#lon_i = numpy.array([-20.0])
#lat_i = numpy.linspace(-90, 90, nlat)  # to regrid to 1/2 degree


# compute cross section interpolation grid
# grid: lat_i, lon_i
# dimensions:  nlat, nlon
if len(lon)*len(lat) == numpy.prod(data2d.shape):
    # lat/lon data.  Find line of data closest to target location
    print("Error: cross section plot not yet supported for structured data")




if var2_read != None:
    dataf = infile.variables[var2_read]

if "ncol_d" in dataf.dimensions or "ncol" in dataf.dimensions:
    print("rank=",dataf.rank,"shape=",dataf.shape,"dims: ",dataf.dimensions)
else:
    print("Error: vertlevels.py only coded for ncol data")

################################################################
# coordinates
################################################################
if "ncol_d" in dataf.dimensions:
    lat  = infile.variables["lat_d"][:]
    lon  = infile.variables["lon_d"][:]
else:
    lat  = infile.variables["lat"][:]
    lon  = infile.variables["lon"][:]

    
################################################################
# time dimension and time values to plot
################################################################
timedim =  "time" in dataf.dimensions
levdim = "lev" in dataf.dimensions or "ilev" in dataf.dimensions
ntimes=1
times=numpy.array([0])
if timedim:
    ntimes = infile.dimensions['time']
    times = infile.variables["time"][:]

if timeindex==None or timeindex==-1:
    t1=ntimes-1      # last timelevel
    t2=ntimes
elif timeindex==-2:
    t1=0             # all timelevels
    t2=ntimes
else:
    t1=timeindex     # user specified index
    t2=timeindex+1

################################################################
# level varaibles 
################################################################
nlev=0

nlev=infile.dimensions["lev"]
lev=infile.variables["lev"][:]
    
hyam=infile.variables['hyam']
hybm=infile.variables['hybm']
hyai=infile.variables['hyai']
hybi=infile.variables['hybi']
print("nlev=",nlev)





#
# read PHIS and interpolate
#
t=t1
# read PHIS, interpolate
if (timedim):
    data2d=dataf[t,...]
    print(t,"time=",times[t])
else:
    data2d=dataf[:]
print("Interpolating unstructured to:",nlat,"x",nlon,"nlat x nlon grid")
zh=numpy.empty([nlev+1,nlat,nlon])   # on interfaces
zh[nlev,:,:] = interp_to_latlon(data2d,lat,lon,lat_i,lon_i)

################################################################
# get correct PS variable, interpolate
################################################################
#ps=infile.variables["ps"]
#ps=infile.variables["DYN_PS"]
#ps=infile.variables["PS"]
#ps_i = interp_to_latlon(ps,lat,lon,lat_i,lon_i)
ps_i = ps0 * numpy.exp( -zh[nlev,:,:]/(Rgas*TREF))


# compute level height
Rgas=287.04
g=9.87
Cp=1005.
TREF=288.
T1=6.5e-3*TREF*Cp/g
T0=TREF-T1
ps0=1000*100
for k in range(nlev-1,-1,-1):
    p = hyam[k]*ps0 + hybm[k]*ps_i
    exner = (p/ps0)**(Rgas/Cp)
    dp =  (hyai[k+1]*ps0 + hybi[k+1]*ps_i) - (hyai[k]*ps0 + hybi[k]*ps_i)
    zh[k,:,:] = zh[k+1,:,:]  + Rgas*dp*(T0 + T1*exner)/(p*g)
    inc = Rgas*dp*(T0 + T1*exner)/(p*g)


for k in range(nlev+1):
    print(k,numpy.amin(zh[k,:,:]),numpy.amax(zh[k,:,:]))
        

wks_type = "pdf"
outname=outname+"."+wks_type
print("MPL output file: ",outname)

# mpl line plot of topography:        


    


