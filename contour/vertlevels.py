#!/usr/bin/env python
#
# vertlevels.py  -i PHIS_file  [-j hyam/ps file]  [-m andes]   PHIS/geos
#
# plot vertical level cross section
#   default: equator
#   -m andes   
#
from __future__ import print_function
import os, sys, numpy
#import Nio
from netCDF4 import Dataset
from plotutils import myargs, extract_level, interp_to_latlon
#from matplotlib import pyplot
import matplotlib.pyplot as plt

inname,inname2,varnames,proj,timeindex,klev,plev,clev,nlatlon_interp,use_ngl,scrip_file,gll_file,se_file \
    = myargs(os.sys.argv)


var1_read=varnames[0]  # phis or geo
units=""
longname=""
if inname2=='':
    inname2=inname

print('PHIS file=',inname)
print('vcoord file=',inname2)
infile = Dataset(file(inname,"r")
infile2 = Dataset(inname2,"r")
outname=inname.split(".nc")[0] + ".topo"


nlat=1
nlon=1024
lat_i = 0
if proj=='andes':
    lat_i = numpy.array([-20.0])
if proj=='himalaya':
    lat_i = numpy.array([30.0])
lon_i = numpy.linspace(-180, 180, nlon,endpoint=False)


#nlat=1024
#nlon=1
#lon_i = numpy.array([-20.0])
#lat_i = numpy.linspace(-90, 90, nlat)  # to regrid to 1/2 degree


################################################################
# read PHI and coordinates
################################################################
if var1_read != None:
    dataf = infile.variables[var1_read][:]
if "ncol_d" in dataf.dimensions:
    lat  = infile.variables["lat_d"][:]
    lon  = infile.variables["lon_d"][:]
else:
    lat  = infile.variables["lat"][:]
    lon  = infile.variables["lon"][:]

print("PHIS: rank=",dataf.rank,"shape=",dataf.shape,"dims: ",dataf.dimensions)
if "ncol_d" in dataf.dimensions or "ncol" in dataf.dimensions:
    idx_lat=None
else:
    idx_lat=( (lat-lat_i)**2).argmin()   # closest point to lat_i
    print("lat_i, lat(idx)",lat_i,lat[idx_lat])
    nlon=len(lon)
    lon_i=lon


    
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

print("reading hybrid coordinate data")
nlev=infile2.dimensions["lev"][:]
lev=infile2.variables["lev"][:]
    
hyam=infile2.variables['hyam'][:]
hybm=infile2.variables['hybm'][:]
hyai=infile2.variables['hyai'][:]
hybi=infile2.variables['hybi'][:]
print("nlev=",nlev)


# constants
Rgas=287.04
g=9.87
Cp=1005.
TREF=288.
T1=6.5e-3*TREF*Cp/g
T0=TREF-T1
ps0=1000*100



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
data2d=data2d/g
print("ZS min,max: ",numpy.amin(data2d),numpy.amax(data2d))
zh=numpy.empty([nlev+1,nlat,nlon])   # on interfaces
if idx_lat==None:
    print("Interpolating unstructured to:",nlat,"x",nlon,"nlat x nlon grid")
    zh[nlev,:,:] = interp_to_latlon(data2d,lat,lon,lat_i,lon_i)
    print("ZS interpoalted min,max: ",numpy.amin(zh[nlev,:,:]),numpy.amax(zh[nlev,:,:]))
else:
    zh[nlev,0,:] = data2d[idx_lat,:]
    print("ZS at idx_lat min,max: ",numpy.amin(zh[nlev,:,:]),numpy.amax(zh[nlev,:,:]))



################################################################
# get correct PS variable, interpolate
################################################################
#ps=infile.variables["ps"]
#ps=infile.variables["DYN_PS"]
#ps=infile.variables["PS"]
#ps_i = interp_to_latlon(ps,lat,lon,lat_i,lon_i)
ps_i = ps0 * numpy.exp( -zh[nlev,:,:]/(Rgas*TREF))


ph=numpy.empty([nlev,nlat,nlon]) 
for k in range(nlev-1,-1,-1):
    p = hyam[k]*ps0 + hybm[k]*ps_i
    exner = (p/ps0)**(Rgas/Cp)
    dp =  (hyai[k+1]*ps0 + hybi[k+1]*ps_i) - (hyai[k]*ps0 + hybi[k]*ps_i)
    zh[k,:,:] = zh[k+1,:,:]  + Rgas*dp*(T0 + T1*exner)/(p*g)
    ph[k,:,:] = p
    inc = Rgas*dp*(T0 + T1*exner)/(p*g)


#for k in range(nlev+1):
#    print(k,numpy.amin(zh[k,:,:]),numpy.amax(zh[k,:,:]))
        

outname="topox.pdf"
print("MPL output file: ",outname)


# mpl line plot of topography:        
fig, axs = plt.subplots()
fig2, axs2 = plt.subplots()

if len(lat_i)==1:
    inc=5
    for k in range(nlev,-1,-inc):
        axs.plot(lon_i,zh[k,0,:],label='k',color='b')
    axs.set(xlabel='longitude', ylabel='height (m)',title='levels')

    inc=5
    for k in range(nlev-1,-1,-inc):
        axs2.plot(lon_i,ph[k,0,:],label='k',color='b')
    axs2.set(xlabel='longitude', ylabel='Pressure (Pa)',title='levels')


if len(lon_i)==1:
    axs.plot(lon_i,zh[nlev,:,0],label='k',color='b')
    axs.set(xlabel='latitude', ylabel='height (m)',title='levels')
    
axs.grid(True)
plt.show()                                                                                               
axs2.grid(True)
plt2.show()                                                                                               
print("writing plot...")
plt.savefig("temp.png")



