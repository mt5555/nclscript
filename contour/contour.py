#!/usr/bin/env python
# Contour over maps of lat/lon data
#
# contour.py -i inputfile  [options]  varname
#
# Should run if Ngl module is not aviable, as long was one is using
# MPL plots, and avoids pressure interpolation in extract_level()
# pressure interpolation requires Ngl.vinth2p (and will have a runtime
# error of Ngl is missing)
#
#from __future__ import print_function
import os, numpy
from netCDF4 import Dataset
from matplotlib import pyplot
from plotutils import mpl_plot, myargs, interp_to_latlon, mpl_streamlines
from nglutils import  ngl_plot, ngl_open, ngl_end, ngl_read_colormap, \
                      extract_level, ngl_vertprofile


inname,inname2,varnames,proj,timeindex,klev,plev,clev,nlatlon_interp,use_ngl, \
scrip_file,gll_file,se_file,contour_opt,coutlines \
    = myargs(os.sys.argv)

var1 = varnames[0]
var1_read=var1
var2_read=None
scale=None
units=""
longname=""
print('file=',inname)
print('contour:',var1,'proj=',proj)
infile = Dataset(inname,"r")
outname=inname.split(".nc")[0] + "."+var1



################################################################
# custom axis for vertical profiles
################################################################
ybnds=[]
xbnds=[]
if var1=="dtheta_dp":
    ybnds=[-.4,.2]
    ybnds=[-.01,.002]
if var1=="dt_dp":
    ybnds=[-.01,.002]
if var1=="Th":
    ybnds=[0,1000]
    xbnds=[150,450]
if var1=="div":
    ybnds=[0,1000]
    xbnds=[-.0008,.0004]

################################################################
# custom contour levels for certain varaiables
################################################################
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
    if var1=="SWCF":
        clev=[-110.,0.,2]
    if var1=="TMQ":
        clev=[0.,50.,.5]
    if var1=="Th":
        clev=[200.,370.,10]
    if var1=="POTT":
        clev=[200.,370.,10]
    if var1=="dtheta_dp":
        clev=[-.1,.1,.005]
    #if var1=="dt_dp":
    #    clev=[-.002,.002,.0001]

else:
    if len(clev)==2:
        clev.append((clev[1]-clev[0])/40)


################################################################
# custom scalings for certain varaiables
################################################################
if var1=="PRECT" or var1=="PRECC" or var1=="PRECL" or var1=="PrecipTotalSurfMassFlux":     
    scale=1000.0*(24*3600)  # convert to mm/day
    units="mm/day"

if var1=="ps":
    scale=1/100.0e0 
    units="hPa"



################################################################
# special processing
################################################################
compute_psdiff=False
if var1=="pnh-ps": 
    print("Special processing:  pnhdiff=pnh-ps") 
    compute_psdiff=True
    var1_read="pnh_i"
    longname="pnhs-ps"

compute_mu=False
if var1=="mu": 
    print("Special processing:  mu = d(PNH)/dp3d") 
    compute_mu=True
    var1_read="DYN_PNH"
    longname="mu"

compute_dz=False
if var1=="DZ3": 
    print("Special processing:  dz = DYN_Z3(k)-DYN_Z3(k+1)") 
    compute_dz=True
    var1_read="DYN_Z3"
    longname="DYN_DZ"

compute_windstress=False
if var1=="surf_mom_flux":
    print("Special processing:  windstress=surf_mom_flux^2") 
    compute_windstress=True


compute_prect=False
if var1=="PRECT" and ~("PRECT" in infile.variables.keys()):
    print("Special processing:  PRECT=PRECC+PRECL") 
    compute_prect=True
    var1_read="PRECL"
    PRECC = infile.variables["PRECC"]
    longname="PRECT"

compute_dtheta_dp=False
if var1=="dtheta_dp":
    print("Special processing:  d(theta)/dp") 
    compute_dtheta_dp=True
    var1_read="Th"
    longname="dtheta/dp"
    units="K/Pa"
    if plev != None:
        print("d(theta)/dp processing only supported on model levels")
        sys.exit(2)


compute_dt_dp=False
if var1=="dt_dp":
    print("Special processing:  dT/dp") 
    compute_dt_dp=True
    var1_read="T"
    longname="dT/dp"
    units="K/Pa"
    if plev != None:
        print("dT/dp processing only supported on model levels")
        sys.exit(2)


compute_pott=False
if var1=="POTT":
    print("Special processing:  POT TEMP") 
    compute_pott=True
    var1_read="T"
    longname="potential temperature"
    units="K"
    if plev != None:
        print("POT TEMP processing only supported on model levels")
        sys.exit(2)


compute_avedx=False
if var1=="ave_dx" or var1=="ave_dx_T":
    print("Special processing:  ave_dx=(min_dx+max_dx)/2") 
    compute_avedx=True
    var1_read="min_dx"
    var2_read="max_dx"
    longname="ave_dx"
    units="km"


compute_streamlines=False
if var1=="stream_uv":
    print("Special processing:  ave_dx=(min_dx+max_dx)/2") 
    compute_streamlines=True
    var1_read="u"
    var2_read="v"
    longname="velocity"
    units="m/s"


    
print("reading dataf name=",var1_read)    
dataf  = infile.variables[var1_read]
data2d_plot2=numpy.array([])
if var2_read != None:
    dataf2 = infile.variables[var2_read]
print("dataf: ",dataf.shape,dataf.dimensions)


title=var1
if longname=="":
    if hasattr(dataf,"long_name"):
        longname=dataf.long_name
    else:
        longname=var1
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
    ntimes = infile.dimensions['time'].size
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

time_label=False
if len(range(t1,t2)) > 1 :
    time_label=True   # add t= label to plots
    print("Adding t= label to plots")

################################################################
# level varaibles
################################################################
nlev=0
nlev_data=0
if levdim:
    nlev=infile.dimensions["lev"].size
    #lev=infile.variables["lev"][:]
    if "lev" in dataf.dimensions:
        nlev_data=nlev
    if "ilev" in dataf.dimensions:
        nlev_data=infile.dimensions["ilev"].size

    if klev==None:
        klev=int(3*nlev/4)

if "hyam" in infile.variables.keys():
    hyam=infile.variables['hyam'][:]
    hybm=infile.variables['hybm'][:]
    hyai=infile.variables['hyai'][:]
    hybi=infile.variables['hybi'][:]
else:
    hyam=None
    hybm=None
    hyai=None
    hybi=None


################################################################
# get correct PS variable
################################################################
ps=numpy.empty([ntimes,1])  # dummy array, wont be used, but needs to be indexed    
have_ps=False
if "ps" in infile.variables.keys():
    ps0=1000*100
    ps=infile.variables["ps"]
    have_ps=True

if "P0" in infile.variables.keys():
    ps0=infile.variables["P0"].getValue()

if "ncol_d" in dataf.dimensions:
    lat  = infile.variables["lat_d"][:]
    lon  = infile.variables["lon_d"][:]
    PSname="DYN_PS"
    if "DYN_PS" in infile.variables.keys():
        ps=infile.variables["DYN_PS"]
        have_ps=True
else:
    lat  = infile.variables["lat"][:]
    lon  = infile.variables["lon"][:]
    PSname="PS"
    if "PS" in infile.variables.keys():
        ps=infile.variables["PS"]
        have_ps=True


if plev != None:
    if not have_ps:
        print("Error: need PS to interpolate to plev=",plev)
        sys.exit(2)
    print("Interpolating to pressure level = ",plev,"using",PSname)


    

#wks_type = "pdf"
wks_type = "png"
    
if use_ngl:
    wks = ngl_open(wks_type,outname)
    if levdim and nlev>1:
        wks_v = ngl_open(wks_type,outname+"_v")
    print("NGL output file: ",outname+"."+wks_type)
    cmap='MPL_viridis'
    #cmap="WhiteBlueGreenYellowRed"
    #cmap="wgne15"
    #cmap="StepSeq25"
    #cmap="BlAqGrYeOrReVi200"
    if compute_avedx:
        cmap = ngl_read_colormap(cmap)
        cmap=numpy.flip(cmap,0)


    if len(clev)==3:
        if clev[1]==-clev[0]:
            cmap='MPL_RdYlBu'
            #cmap="BlueYellowRed"

    if var1=="ps" and len(clev)==3:
        # custom colormap with center at 1000mb
        cmap = ngl_read_colormap("MPL_RdYlBu")
        n=cmap.shape
        n=n[0]-1
        nint=(clev[1]-clev[0])/clev[2]
        print("cmap=",n,cmap.shape)
        if (1000-clev[0]) > (clev[1]-1000):
            x1=0
            x2=int(n/nint +    ( (clev[1]-clev[0]) / (1000-clev[0]) )*n/2 )
        else:
            x1=int(n/nint +  n-  ( (clev[1]-clev[0]) / (clev[1]-1000) )*n/2 )
            x2=n
        print("using symmetric contour map centered at 1000mb:",cmap.shape,x1,x2)
        cmap=cmap[x1:x2,:]


else:
    outname=outname+"."+wks_type
    print("MPL output file: ",outname)
    #cmap='nipy_spectral'
    #cmap='viridis'
    cmap='plasma'
    if len(clev)==3:
        if clev[1]==-clev[0]:
            cmap='Spectral'     # good diverging colormap
            #cmap='RdYlBu'     # good diverging colormap
    if var1=="ps" and len(clev)==3:
        cmap='RdYlBu'

for t in range(t1,t2):

    #
    # 2D maps
    #
    if (timedim and levdim):
        dimname=dataf.dimensions
        if "lev" in dimname:
            kidx=dimname.index("lev")
        if "ilev" in dimname:
            kidx=dimname.index("ilev")
        print(t+1,"time=",times[t],"k=",klev+1,"/",nlev_data,"plev=",plev," level index=",kidx)
        data2d=extract_level(dataf[t,...],klev,plev,ps[t,...],hyam,hybm,kidx-1)
        print("data2d shape=",data2d.shape)

        if compute_streamlines:
            data2d_2=extract_level(dataf2[t,...],klev,plev,ps[t,...],hyam,hybm,kidx-1)            

        if compute_dtheta_dp or compute_dt_dp:
            data2dm1=extract_level(dataf[t,...],klev-1,plev,ps[t,...],hyam,hybm)
            # midpoint data
            dp =  (hyam[klev]*ps0 + hybm[klev]*ps[t,...]) -  \
                (hyam[klev-1]*ps0 + hybm[klev-1]*ps[t,...]) 
            data2d = (data2d - data2dm1)/dp

        if compute_dz:
            data2dm1=extract_level(dataf[t,...],klev+1,plev,ps[t,...],hyam,hybm)
            data2d = (data2d-data2dm1)

        if compute_pott:
            p =  hyam[klev]*ps0 + hybm[klev]*ps[t,...]
            exner=(p/ps0)**.2856  # kappa
            data2d=data2d/exner

        if compute_psdiff:
            data2d=data2d - ps[t,...]

        if compute_mu:            
            # interface data:
            dp =  (hyai[klev+1]*ps0 + hybi[klev+1]*ps[t,...]) -  \
                (hyai[klev]*ps0 + hybi[klev]*ps[t,...]) 
            data2dm1=extract_level(dataf[t,...],klev+1,plev,ps[t,...],hyam,hybm)
            data2d=(data2dm1-data2d)/dp
            


    elif (timedim):
        print(t+1,"time=",times[t]," time & space dimensions only")
        dimname=dataf.dimensions
        #data2d=numpy.squeeze(dataf[t,...],0)
        data2d=dataf[t,...]
        if compute_prect:
            data2d=data2d + PRECC[t,...]
            longname="PRECT"
        if compute_windstress:
            print("fix this code - should this be u^2 + v^2?")
            sys.exit(1)
            temp=data2d[:,0]**2 + data2d[:,0]**2
            data2d=temp
            longname="wind stress"
    elif (levdim):
        print("No time dimension.  k=",klev,"/",nlev_data,"plev=",plev)
        data2d=extract_level(dataf[...],klev,plev,ps[...],hyam,hybm)
    else:
        print("No time or level dimension, assuming space dimension only")
        data2d=dataf[:]
        if compute_avedx:
            data2d_2=dataf2[:]
            data2d=numpy.sqrt(data2d*data2d_2)
            if var1=="ave_dx_T":
                dataf_T  = infile.variables["T"]
                print("reading T to add to avedx plot",dataf_T.shape,ps.shape)
                data2d_plot2 = extract_level(dataf_T[6,...],klev,[750.0],ps[6,...],hyam,hybm)
            else:
                # plot contour lines of ps on top of avedx:
                data2d_plot2=ps[t,...]

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

    if time_label:
        title=longname
        longname_t="t="+str(times[t])
    else:
        longname_t=longname

    # should we interpolate?
    if len(lon)*len(lat) != numpy.prod(data2d.shape) and nlatlon_interp:
        nlat=nlatlon_interp[0]
        nlon=nlatlon_interp[1]
        # dont add endpoint, to be consistent with lat x lon history files
        lon_i = numpy.linspace(0, 360, nlon,endpoint=False)  
        if nlat % 2 == 0:
            dlat2=90./nlat
            lat_i = numpy.linspace(-90+dlat2, 90-dlat2, nlat)  
            print("Interpolating unstructured to:",nlat,"x",nlon,"nlat x nlon uni grid")
        else:
            lat_i = numpy.linspace(-90, 90, nlat)  # to regrid to 1/2 degree
            print("Interpolating unstructured to:",nlat,"x",nlon,"nlat x nlon cap grid")
            
        data_i=interp_to_latlon(data2d,lat,lon,lat_i,lon_i)
        if use_ngl:
            ngl_plot(wks,data_i,lon_i,lat_i,title,longname_t,units,
                     proj,clev,cmap,scrip_file,se_file,contour_opt,coutlines)
        else:
            if compute_streamlines:
                data_i_2=interp_to_latlon(data2d_2,lat,lon,lat_i,lon_i)
                mpl_streamlines(data_i,data_i_2,lon_i,lat_i,title,longname,units,
                         proj,clev,cmap)
            else:
                print("calling mpl_plot")
                mpl_plot(data_i,lon_i,lat_i,title,longname,units,
                         proj,clev,cmap,scrip_file,gll_file,contour_opt,coutlines)
            pyplot.savefig(outname,dpi=300,orientation="portrait")
            pyplot.close()
    elif use_ngl:
        ngl_plot(wks,data2d,lon,lat,title,longname_t,units,
                 proj,clev,cmap,scrip_file,se_file,contour_opt,coutlines,data2d_plot2)
    else:
        if compute_streamlines:
            mpl_streamlines(data2d,data2d_2,lon,lat,title,longname,units,
                            proj,clev,cmap)
        else:
            mpl_plot(data2d,lon,lat,title,longname,units,
                     proj,clev,cmap,scrip_file,gll_file,contour_opt,coutlines)
        #pyplot.show()
        pyplot.savefig(outname,dpi=300,orientation="portrait")
        pyplot.close()


    #
    # 1D vertical profile at a specified point
    #
    if levdim and nlev_data>1 and have_ps and use_ngl:
        # data2d[min_i1] will give value for either 1D or 2D data
        min_i1=[183,115]
        max_i1=[183,112]
        latmin  = lat[min_i1[0]]
        lonmin  = lon[min_i1[1]]
        latmax  = lat[max_i1[0]]
        lonmax  = lon[max_i1[1]]
        print("vert profile1 at ",min_i1,"lon=",lonmin,'lat=',latmin)
        print("vert profile2 at ",max_i1,"lon=",lonmax,'lat=',latmax)

        ncols=len((min_i1,max_i1))
        # add ref profile to plot?
        #if var1=="Th" or var1=="POTT":
        #    ncols=ncols*3
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

            if compute_dtheta_dp or compute_dt_dp:
                dp = p_i[1:nlev+1]-p_i[0:nlev]
                # compute d/dp
                coldata_i=p_i
                coldata_i[0]=coldata[0]
                coldata_i[nlev]=coldata[nlev-1]
                coldata_i[1:nlev] = (coldata[0:nlev-1]+coldata[1:nlev])/2
                coldata = (coldata_i[1:nlev+1]-coldata_i[0:nlev])/dp

            if compute_pott:
                exner=(p/ps0)**.2856  # kappa
                coldata=coldata/exner
                
                
            coldata_all[:,icols]=coldata
            if (nlev_data==nlev):
                p_all[:,icols]=p
            else:
                p_all[:,icols]=p_i
            icols=icols+1
            if False and (var1=="Th" or var1=="POTT"):
                print("Adding Th reference profiles")
                exner=(p/ps0)**.2856  # kappa
                coldata_all[:,icols]=97/exner + 191   # ECMWF
                p_all[:,icols]=p
                icols=icols+1
                coldata_all[:,icols]=300*exner**(0.67-1)  # Shar
                p_all[:,icols]=p
                icols=icols+1

        kstart=int(nlev/3)
        ngl_vertprofile(wks_v,coldata_all[kstart:,...],p_all[kstart:,...],xbnds,ybnds,longname,units,times[t])

    


if use_ngl:
    del wks
    if "wks_v" in locals():
        del wks_v
    ngl_end()
