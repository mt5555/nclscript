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
from vertprofile import ngl_vertprofile

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

outname=inname.split(".nc")[0] + "."+var1

vrange=[]
if var1=="dtheta_dp":
    vrange=[-.4,.2]

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
    if var1=="Th":
        clev=[200.,370.,10]
    if var1=="dtheta_dp":
        clev=[-.1,.1,.005]

else:
    if len(clev)==2:
        clev.append((clev[1]-clev[0])/40)

infile = Nio.open_file(inname,"r")

################################################################
# special processing
################################################################
PRECT_from_LC=False
if var1=="PRECT" and ~("PRECT" in infile.variables.keys()):
    print("Special processing:  PRECT=PRECC+PRECL") 
    PRECT_from_LC=True
    var1_read="PRECL"
    PRECC = infile.variables["PRECC"]
    longname="PRECT"

if var1=="PRECT" or var1=="PRECC" or var1=="PRECL":     
    scale=1000.0*(24*3600)  # convert to mm/day
    units="mm/day"

compute_dtheta_dp=False
if var1=="dtheta_dp":
    print("Special processing:  d(theta)/dp") 
    compute_dtheta_dp=True
    var1_read="Th"
    longname="dtheta/dp"
    units="K/Pa"



    
dataf  = infile.variables[var1_read]
print("rank=",dataf.rank,"shape=",dataf.shape,"dims: ",dataf.dimensions)




title=var1
if longname=="" and hasattr(dataf,"long_name"):
    longname=dataf.long_name
    title=""
if units=="" and hasattr(dataf,"units"):
    units=dataf.units
        
# does array have a time dimension?
timedim =  "time" in dataf.dimensions
levdim = "lev" in dataf.dimensions
ntimes=1
times=numpy.array([0])
if timedim:
    ntimes = infile.dimensions['time']
    times = infile.variables["time"][:]


#
# get correct PS variable
#
ps=numpy.empty([ntimes,1])  # dummy array, wont be used, but needs to be indexed    
if "ps" in infile.variables.keys():
    ps0=1000*100
    ps=infile.variables["ps"]

if "ncol_d" in dataf.dimensions:
    lat  = infile.variables["lat_d"][:]
    lon  = infile.variables["lon_d"][:]
    ps0=infile.variables["P0"].get_value()
    PSname="DYN_PS"
    ps=infile.variables["DYN_PS"]

else:
    lat  = infile.variables["lat"][:]
    lon  = infile.variables["lon"][:]
    ps0=infile.variables["P0"].get_value()
    PSname="PS"
    ps=infile.variables["PS"]

if plev != None:
    print("Interpolating to pressure level = ",plev,"using",PSname)
    klev=-1  # flag indicating interpolation


nlev=0
if "lev" in dataf.dimensions:
    nlev=infile.dimensions["lev"]
    lev=infile.variables["lev"][:]
    if klev==None:
        klev=int(3*nlev/4)
    hyam=infile.variables['hyam']
    hybm=infile.variables['hybm']
    hyai=infile.variables['hyai']
    hybi=infile.variables['hybi']


        


# default plot all times:

if timeindex==None or timeindex==-1:
    t1=ntimes-1      # last timelevel
    t2=ntimes
elif timeindex==-2:
    t1=0             # all timelevels
    t2=ntimes
else:
    t1=timeindex     # user specified index
    t2=timeindex+1



    

wks_type = "pdf"
    
if use_ngl:
    wks = Ngl.open_wks(wks_type,outname)
    if levdim and nlev>0:
        wks_v = Ngl.open_wks(wks_type,outname+"_v")
    print("output file: ",outname+"."+wks_type)
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
    outname=outname+"."+wks_type
    print("output file: ",outname)
    #cmap='nipy_spectral'
    #cmap='viridis'
    cmap='plasma'
    if len(clev)==3:
        if clev[1]==-clev[0]:
            cmap='Spectral'     # good diverging colormap
            #cmap='RdYlBu'     # good diverging colormap
    

for t in range(t1,t2):

    #
    # 2D maps
    #
    if (timedim and levdim):
        print(t+1,"time=",times[t],"k=",klev+1,"/",nlev,"plev=",plev)
        data2d=extract_level(dataf[t,...],klev,plev,ps[t,...],hyam,hybm)
        if compute_dtheta_dp:
            data2dm1=extract_level(dataf[t,...],klev-1,plev,ps[t,...],hyam,hybm)
            dp =  (hyam[klev]*ps0 + hybm[klev]*ps[t,...]) -  \
                (hyam[klev-1]*ps0 + hybm[klev-1]*ps[t,...]) 
            data2d = (data2d - data2dm1)/dp

    elif (timedim):
        data2d=dataf[t,...]
        if PRECT_from_LC:
            data3d=data2d + PRECC[t,...]
            longname="PRECT"
        print(t,"time=",times[t])
    elif (levdim):
        print("k=",klev,"/",nlev,"plev=",plev)
        data2d=extract_level(dataf[...],klev,plev,ps[...],hyam,hybm)
    else:
        data2d=dataf[:]

    if scale:
        data2d=data2d*scale


    data2dmin=numpy.amin(data2d)
    data2dmax=numpy.amax(data2d)
    min_i = numpy.where(data2d == data2dmin)
    max_i = numpy.where(data2d == data2dmax)
    # take the first occurance:
    if len(data2d.shape)==2:
        min_i1=[min_i[0][0],min_i[1][0]]
        max_i1=[max_i[0][0],max_i[1][0]]
        latmin  = lat[min_i1[0]]
        lonmin  = lon[min_i1[1]]
        latmax  = lat[max_i1[0]]
        lonmax  = lon[max_i1[1]]
    else:
        min_i1=[min_i[0][0]]
        max_i1=[max_i[0][0]]
        latmin  = lat[min_i1]
        lonmin  = lon[min_i1]
        latmax  = lat[max_i1]
        lonmax  = lon[max_i1]

    print("scaled data2d min=",data2dmin,"at",min_i1,"lon=",lonmin,'lat=',latmin)
    print("scaled data2d max=",data2dmax,"at",max_i1,"lon=",lonmax,'lat=',latmax)


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


    #
    # 1D vertical profile at a specified point
    #
    # data2d[min_i1] will give value for either 1D or 2D data
    if levdim and nlev>0:
        ncols=len((min_i1,max_i1))
        if var1=="Th":
            ncols=ncols*2  # add ref profile
        coldata_all=numpy.zeros([nlev,ncols])
        p_all=numpy.zeros([nlev,ncols])
        icols=0
        for idx in (min_i1,max_i1):
            if (timedim):
                if len(idx)==2:
                    coldata=dataf[t,:,idx[0],idx[1]]
                else:
                    coldata=dataf[t,:,idx[0]]
            elif (levdim):
                if len(idx)==2:
                    coldata=dataf[:,idx[0],idx[1]]
                else:
                    coldata=dataf[:,idx[0]]
            if scale:
                coldata=coldata*scale

            if len(idx)==2:
                p =  hyam[:]*ps0 + hybm[:]*ps[t,idx[0],idx[1]]
                p_i= hyai[:]*ps0 + hybi[:]*ps[t,idx[0],idx[1]]
            else:
                p =  hyam[:]*ps0 + hybm[:]*ps[t,idx[0]]
                p_i= hyai[:]*ps0 + hybi[:]*ps[t,idx[0]]

            if compute_dtheta_dp:
                dp = p_i[1:nlev+1]-p_i[0:nlev]
                # compute d/dp
                coldata_i=p_i
                coldata_i[0]=coldata[0]
                coldata_i[nlev]=coldata[nlev-1]
                coldata_i[1:nlev] = (coldata[0:nlev-1]+coldata[1:nlev])/2
                coldata = (coldata_i[1:nlev+1]-coldata_i[0:nlev])/dp
                

            coldata_all[:,icols]=coldata
            p_all[:,icols]=p
            icols=icols+1
            if var1=="Th":
                exner=(p/ps0)**.2856  # kappa
                coldata_all[:,icols]=97/exner + 191
                p_all[:,icols]=p
                icols=icols+1

        kstart=int(nlev/3)
        ngl_vertprofile(wks_v,coldata_all[kstart:,...],p_all[kstart:,...],vrange,longname,units)

    


if use_ngl:
    del wks
    if "wks_v" in locals():
        del wks_v
    Ngl.end()
