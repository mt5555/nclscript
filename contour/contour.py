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
from plotutils import mpl_plot, ngl_plot, myargs, extract_level, interp_to_latlon
from matplotlib import pyplot

inname,varnames,proj,timeindex,klev,plev,clev,nlatlon_interp,use_ngl,scrip_file,gll_file \
    = myargs(os.sys.argv)

var1 = varnames[0]
var1_read=var1
var2_read=None
scale=None
units=""
longname=""

print('file=',inname)
print('contour:',var1,'proj=',proj)


if clev==None:
    clev=[40]   # 40 levels, no range specified

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
else:
    if len(clev)==2:
        clev.append((clev[1]-clev[0])/40)

infile = Nio.open_file(inname,"r")


# special processing
PRECT_from_LC=False
if var1=="PRECT" and ~("PRECT" in infile.variables.keys()):
    print("Special processing:  PRECT=PRECC+PRECL") 
    PRECT_from_LC=True
    var1_read="PRECL"
    var2_read="PRECC"
    longname="PRECT"

if var1=="PRECT" or var1=="PRECC" or var1=="PRECL":     
    scale=1000.0*(24*3600)  # convert to mm/day
    units="mm/day"


    
dataf  = infile.variables[var1_read]
print("rank=",dataf.rank,"shape=",dataf.shape,"dims: ",dataf.dimensions)

if var2_read:
    data2 = infile.variables[var2_read]


title=var1
if longname=="" and hasattr(dataf,"long_name"):
    longname=dataf.long_name
    title=""
if units=="" and hasattr(dataf,"units"):
    units=dataf.units
        
    
if "ncol_d" in dataf.dimensions:
    lat  = infile.variables["lat_d"][:]
    lon  = infile.variables["lon_d"][:]
    PSname = "DYN_PS"
else:
    lat  = infile.variables["lat"][:]
    lon  = infile.variables["lon"][:]
    PSname = "PS"


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
    print("Interpolating to pressure level = ",plev,"using",PSname)
    klev=-1  # flag indicating interpolation
    PS=infile.variables[PSname]
    hyam=infile.variables['hyam']
    hybm=infile.variables['hybm']

    

wks_type = "pdf"
    
if use_ngl:
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
    outname=var1+"."+wks_type
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
        if PRECT_from_LC:
            data3d=data2d + data2[t,...]
            longname="PRECT"
        print(t,"time=",times[t])
    elif (levdim):
        print("k=",klev,"/",nlev,"plev=",plev)
        data2d=extract_level(dataf[...],klev,plev,PS[...],hyam,hybm)
    else:
        data2d=dataf.values

    if scale:
        data2d=data2d*scale

    print("scaled data2d min/max=",numpy.amin(data2d),numpy.amax(data2d))

    # should we interpolate?
    if len(lon)*len(lat) != numpy.prod(data2d.shape) and nlatlon_interp:
        nlat=nlatlon_interp[0]
        nlon=nlatlon_interp[1]
        lon_i = numpy.linspace(0, 360, nlon,endpoint=False)  
        if nlat % 2 == 0:
            dlat2=90./nlat
            lat_i = numpy.linspace(-90+dlat2, 90-dlat2, nlat)  
            print("Interpolating unstructured to: ",nlat,"x",nlon,"uni grid")
        else:
            lat_i = numpy.linspace(-90, 90, nlat)  # to regrid to 1/2 degree
            print("Interpolating unstructured to: ",nlat,"x",nlon,"cap grid")
            
        data_i=interp_to_latlon(data2d,lat,lon,lat_i,lon_i)
        ngl_plot(wks,data_i,lon_i,lat_i,title,longname,units,
                 proj,clev,cmap,scrip_file)
    elif use_ngl:
        ngl_plot(wks,data2d,lon,lat,title,longname,units,
                 proj,clev,cmap,scrip_file)
    else:
        mpl_plot(data2d,lon,lat,title,longname,units,
                 proj,clev,cmap,gll_file)
        pyplot.savefig(outname,dpi=300,orientation="portrait")
        #pyplot.show()

if use_ngl:
    del wks
    Ngl.end()
