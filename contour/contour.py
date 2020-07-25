#!/usr/bin/env python
# Contour over maps of lat/lon data
# python contour_latlon.py
#
#
# contourngl -i inputfile  [options]  varname
#
from __future__ import print_function
import os, numpy
import Nio, Ngl
from plotutils import mpl_plot, ngl_plot, myargs, extract_level
from matplotlib import pyplot


inname,varnames,proj,timeindex,klev,plev,clev,use_ngl,scrip_file,gll_file \
    = myargs(os.sys.argv)

var1 = varnames[0]
print('file=',inname)
print('contour:',var1,'proj=',proj)

if clev==None:
    clev=[50]   # 50 levels, no range specified

    if var1=="TBOT":
        clev=[160.,300.,5.]
    if var1=="PRECT":
        clev=[0.,10.,0.25]
    if var1=="PSL":
        clev=[96000.,105000.,200.]
    if var1=="PS":
        clev=[52000.,105000.,2000.]
    if var1=="PHIS":
        clev=[-500.,55000.,1000.]
    if var1=="TMQ":
        clev=[0.,50.,.5]

        

infile = Nio.open_file(inname,"r")
dataf  = infile.variables[var1]
print("rank=",dataf.rank,"shape=",dataf.shape,"dims: ",dataf.dimensions)

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


if "lev" in dataf.dimensions:
    nlev=infile.dimensions["lev"]
    if klev==None:
        klev=int(3*nlev/4)

        
# does array have a time dimension?
ntimes=1
times=numpy.array([0])
if "time" in dataf.dimensions:
    ntimes = infile.dimensions['time']
    times = infile.variables["time"][:]

# default plot all times:

if timeindex==None:
    t1=ntimes-1      # last timelevel
    t2=ntimes
elif timeindex==-1:
    t1=0             # all timelevels
    t2=ntimes
else:
    t1=timeindex     # user specified index
    t2=timeindex+1


if plev == None:
    PS=numpy.empty([ntimes,1])  # dummy array, wont be used, but needs to be indexed
    hyam=None
    hybm=None
else:
    print("Interpolating to pressure level = ",plev)
    klev=-1  # flag indicating interpolation
    PS=infile.variables['DYN_PS']
    hyam=infile.variables['hyam']
    hybm=infile.variables['hyam']
    
            

    
if use_ngl:
    wks_type = "png"
    wks = Ngl.open_wks(wks_type,var1)
    cmap='MPL_viridis'
    #cmap="WhiteBlueGreenYellowRed"
    #cmap="wgne15"
    #cmap="StepSeq25"
    #cmap="BlAqGrYeOrReVi200"
    if len(clev)==3:
        if clev[1]==-clev[0]:
            cmap='MPL_RdYlBu'
            #cmap="BlueYellowRed"


else:
    outname=var1+".png"
    #cmap='nipy_spectral'
    #cmap='viridis'
    cmap='plasma'
    if len(clev)==3:
        if clev[1]==-clev[0]:
            cmap='Spectral'     # good diverging colormap
            #cmap='RdYlBu'     # good diverging colormap
    

for t in range(t1,t2):
    timedim =  "time" in dataf.dimensions
    levdim = "lev" in dataf.dimensions

    if (timedim and levdim):
        print(t,"time=",times[t],"k=",klev,"/",nlev,"plev=",plev)
        data2d=extract_level(dataf[t,...],klev,plev,PS[t,...],hyam,hybm)
    elif (timedim):
        data2d=dataf[t,...]
        print(t,"time=",times[t])
    elif (levdim):
        print("k=",klev,"/",nlev,"plev=",plev)
        data2d=extract_level(dataf[...],klev,plev,PS[...],hyam,hybm)
    else:
        data2d=dataf.values

    if use_ngl:
        ngl_plot(wks,wks_type,data2d,lon,lat,title,longname,units,
                 proj,clev,cmap,scrip_file)
    else:
        mpl_plot(data2d,lon,lat,title,longname,units,
                 proj,clev,cmap,gll_file)
        pyplot.savefig(outname,dpi=300,orientation="portrait")
        #pyplot.show()

if use_ngl:
    del wks
    Ngl.end()
