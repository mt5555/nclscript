#!/usr/bin/env python
# Contour over maps of lat/lon data
# python contour_latlon.py
#
#
# contour.py -i inputfile  [options]  varname
#
from __future__ import print_function
import os, numpy
import Nio, Ngl
from plotutils import mpl_plot, ngl_plot, myargs, extract_level, interp_to_latlon
from matplotlib import pyplot
from vertprofile import ngl_vertprofile

inname,inname2,varnames,proj,timeindex,klev,plev,clev,nlatlon_interp,use_ngl,scrip_file,gll_file \
    = myargs(os.sys.argv)

var1 = "T"
var1_read=var1
var2_read=None
scale=None
units=""
longname=""
print('file=',inname)
print('contour:',var1,'proj=',proj)
infile = Nio.open_file(inname,"r")
outname=inname.split(".nc")[0] + "."+var1



################################################################
# custom contour levels for certain varaiables
################################################################
vrange=[]

if clev==None:
    clev=[30000.,60000.,1000]
else:
    if len(clev)==2:
        clev.append((clev[1]-clev[0])/40)


dataf  = infile.variables[var1_read]
if var2_read != None:
    dataf2 = infile.variables[var2_read]
print("rank=",dataf.rank,"shape=",dataf.shape,"dims: ",dataf.dimensions)

# todo: read in PHIS, add to Z
title="Z-ZS"
longname="Z-ZS"
units="m"

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
nlev_data=0
if levdim:
    nlev=infile.dimensions["lev"]
    lev=infile.variables["lev"][:]
    if "lev" in dataf.dimensions:
        nlev_data=nlev
    if "ilev" in dataf.dimensions:
        nlev_data=infile.dimensions["ilev"]

    if klev==None:
        klev=int(3*nlev/4)
    hyam=infile.variables['hyam']
    hybm=infile.variables['hybm']
    hyai=infile.variables['hyai']
    hybi=infile.variables['hybi']



################################################################
# get correct PS variable
################################################################
ps=numpy.empty([ntimes,1])  # dummy array, wont be used, but needs to be indexed    
if "ps" in infile.variables.keys():
    ps0=1000*100
    ps=infile.variables["ps"]

if "P0" in infile.variables.keys():
    ps0=infile.variables["P0"].get_value()

if "ncol_d" in dataf.dimensions:
    lat  = infile.variables["lat_d"][:]
    lon  = infile.variables["lon_d"][:]
    PSname="DYN_PS"
    ps=infile.variables["DYN_PS"]
else:
    lat  = infile.variables["lat"][:]
    lon  = infile.variables["lon"][:]
    PSname="PS"
    if "PS" in infile.variables.keys():
        ps=infile.variables["PS"]


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
    
Rgas=287.04
g=9.87
for t in range(t1,t2):

    # sum over levels
    zh=0*ps[t,...]
    for klev in reversed(range(nlev_data)):
        print(t+1,"time=",times[t],"k=",klev+1,"/",nlev_data,"plev=",plev)
        T2d=extract_level(dataf[t,...],klev,plev,ps[t,...],hyam,hybm)

        p = hyam[klev]*ps0 + hybm[klev]*ps[t,...]
        dp =  (hyai[klev+1]*ps0 + hybi[klev+1]*ps[t,...]) -  \
            (hyai[klev]*ps0 + hybi[klev]*ps[t,...])
        zh = zh+ Rgas*dp*T2d/(p*g)
        # need to add PHIS before this data makes sense, dont print for now:
        #print("k_i=",klev+1,"pmin,pmax",numpy.amin(p)/100,numpy.amax(p)/100,"z min,max=",numpy.amin(zh),numpy.amax(zh))
        # on the equator, over the pacific. lat=0, lon=180
        idx=(lat**2+(lon-180)**2).argmin()
        print("k_i=",klev,"At Pacfific Equator p=",p[idx]/100,"mb z=",zh[idx],"m")

    data2d=zh

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
            print("Interpolating unstructured to:",nlat,"x",nlon,"nlat x nlon uni grid")
        else:
            lat_i = numpy.linspace(-90, 90, nlat)  # to regrid to 1/2 degree
            print("Interpolating unstructured to:",nlat,"x",nlon,"nlat x nlon cap grid")
            
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
    if levdim and nlev_data>0:
        ncols=len((min_i1,max_i1))
        coldata_all=numpy.zeros([nlev_data,ncols])
        p_all=numpy.zeros([nlev_data,ncols])
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

                
            coldata_all[:,icols]=coldata
            if (nlev_data==nlev):
                p_all[:,icols]=p
            else:
                p_all[:,icols]=p_i
            icols=icols+1


        kstart=int(nlev/3)
        ngl_vertprofile(wks_v,coldata_all[kstart:,...],p_all[kstart:,...],vrange,longname,units)

    


if use_ngl:
    del wks
    if "wks_v" in locals():
        del wks_v
    Ngl.end()
