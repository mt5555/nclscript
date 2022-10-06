#!/usr/bin/env python
#
# diffps.py  -i input1  -j input2  [options]  varname
#
# plot map of input2:varname - input1:varname
# also compute l2 norm
#
#
from __future__ import print_function
import os, numpy
import Nio, Ngl
from plotutils import mpl_plot, ngl_plot, myargs, extract_level, interp_to_latlon
from matplotlib import pyplot

inname,inname2,varnames,proj,timeindex,klev,plev,clev,nlatlon_interp,use_ngl,scrip_file,gll_file \
    = myargs(os.sys.argv)

var1 = varnames[0]
#var1_read=var1
scale=None
units=""
longname=""
print('file=',inname)
print('contour:',var1,'proj=',proj)
infile = Nio.open_file(inname,"r")
infile2 = Nio.open_file(inname2,"r")
outname=inname.split(".nc")[0] + "."+var1+".diff"



################################################################
# contour levels 
################################################################
if clev==None:
    clev=[-100.,100.,5.]
else:
    if len(clev)==2:
        clev.append((clev[1]-clev[0])/40)


dataf  = infile.variables[var1]
print("file1: rank=",dataf.rank,"shape=",dataf.shape,"dims: ",dataf.dimensions)
dataf2  = infile2.variables[var1]
print("file2: rank=",dataf2.rank,"shape=",dataf2.shape,"dims: ",dataf2.dimensions)

title=var1
if longname=="" and hasattr(dataf,"long_name"):
    longname=dataf.long_name
    title=""
if units=="" and hasattr(dataf,"units"):
    units=dataf.units

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

if plev != None:
    print("Interpolating to pressure level = ",plev,"using",PSname)



    

wks_type = "pdf"
    
if use_ngl:
    wks = Ngl.open_wks(wks_type,outname)
    print("output file: ",outname+"."+wks_type)
    cmap='MPL_viridis'
    if len(clev)==3:
        if clev[1]==-clev[0]:
            cmap='MPL_RdYlBu'
            #cmap="BlueYellowRed"


else:
    outname=outname+"."+wks_type
    print("output file: ",outname)
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
        print(t+1,"time=",times[t],"k=",klev+1,"/",nlev_data,"plev=",plev)
        data2d_1=extract_level(dataf[t,...],klev,plev,ps[t,...],hyam,hybm)
        data2d_2=extract_level(dataf2[t,...],klev,plev,ps[t,...],hyam,hybm)
        

    elif (timedim):
        data2d_1=dataf[t,...]
        print(t,"time=",times[t])
    elif (levdim):
        print("k=",klev,"/",nlev_data,"plev=",plev)
        data2d_1=extract_level(dataf[...],klev,plev,ps[...],hyam,hybm)
    else:
        data2d_1=dataf[:]

    if scale:
        data2d_1=data2d*scale
        data2d_2=data2d_2*scale

    # plot the diff
    data2d=data2d_2-data2d_1

    # l2 norm      
    if len(lon)*len(lat) == numpy.prod(data2d.shape) and  "gw" in infile.variables.keys():
        gw  = infile.variables["gw"][:]
        tmp1 = (  numpy.sum(gw.dot(data2d_1*data2d_1))/sum(gw)/len(lon)  )**0.5
        tmp2 = (  numpy.sum(gw.dot(data2d_2*data2d_2))/sum(gw)/len(lon)  )**0.5
        tmpd = (  numpy.sum(gw.dot(data2d*data2d))/sum(gw)/len(lon)  )**0.5

        print("l2 norm data1=",tmp1)
        print("l2 norm data2=",tmp2)
        print("l2 norm diff abs,rel=",tmpd,tmpd/tmp1)
        
        tmp1=numpy.amax(abs(data2d_1))
        tmpd=numpy.amax(abs(data2d))
        print("max norm diff abs,rel=",tmpd,tmpd/tmp1)        
    else:
        print("unstructured data - add code here to compute l2 norm")


  






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



if use_ngl:
    del wks
    Ngl.end()
