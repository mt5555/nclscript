import numpy, os, sys, getopt

# needed for mpl_plot
from cartopy import crs
from cartopy.util import add_cyclic_point
from matplotlib import pyplot
from scipy.interpolate import griddata

# needed for ngl_plot
import Ngl
from netCDF4 import Dataset

#
# Input:
#   2D field:
#     data(ncol), lon(ncol), lat(ncol)       unstructured
#     data(nlat,nlon), lon(nlon), lat(nlat)  structured grid data
#   title
#   projection      string denoting projection and region
#                   ex: "latlon","US1","andes","oro", etc...
#
#   [nlevels] or [vmin,vmax,nlevels]     contour bounds and number of contours 
#   cmap            colormap 
#
#   scrip_file      option for NGL plots, user specified dual grid
#   gll_file        option for matplotlib, user specified subcells for triangulation
#   se_file         option to NGL plots, to draw spectral element mesh on top of output
# 
#
# can handle many types of plots:
#
# 1. latlon data:  data(nlat,nlon)  NGL and MPL.   best results.
# 2. unstructured  dual grid:  NGL only
#                  2a: construct dual grid via triangulation.  UGLY
#                  2b: specify dual grid (-s scripfile.nc).  BEST FOR PG2 data
#                      can shade each PG2 cell
#                      but does it still have white borders around each cell?
# 3. unstructured: vertex data:  MPL only
#                  3a: construct triangulation in plot coordinates
#                      gouraud shading of vertex data.  BEST FOR GLL DATA
#                      how does it work for PG2 data?
#                  3b: construct triangulation from GLL subcells (-g subcell.nc)
#                      gouraud shading of vertex data (not as good as 3a)
#                  3c: solid fill each scrip cell (-s scripfile.nc)
#                        PolyCollection approach from Ben Hillman
#

def myargs(argv):
    inputfile = ''
    inputfile2 = ''
    scripfile = ''
    gllfile = ''
    se_file = ''
    varname = ''
    contour_opt = ''
    use_ngl = True
    timeindex = None
    levindex = None
    pressurelev = None
    nlatlon_interp=None
    clev = None
    name = argv[0]
    projection='latlon'
    try:
        opts, args = getopt.getopt(argv[1:],"i:j:s:g:t:k:p:y:c:m:r:e:f:o:")
    except getopt.GetoptError:
        print (name,' -i inputfile [options] varname')
        print (name,' -c nlevels  number of contour levels (ignored in MPL)')
        print (name,' -c cmin,cmax       contour level min,max with 40 levels')
        print (name,' -c cmin,cmax,cinc  contour level min,max,spacing')
        print (name,' -c cmin,cmax,cinc,logbase  log levels, min,max,inc (in log space),base (2 or 10)')
        print (name,' -e Exodus.g file for plotting spectral elements')
        print (name,' -f contour fill: area,raster,la (lines+area), lo (lines only)')
        print (name,' -g gll_subcell_file')
        print (name,' -j inputfile2')
        print (name,' -k levindex starting at 1  [default: 3*nlev/4]')
        print (name,' -m map projeciton:  latlon,US1,oro,andes,hamalaya,etc...')
        print (name,' -o 0,1  continential outlines 1=on/0=off')
        print (name,' -p pressure(mb)  interpolate to pressure level')
        print (name,' -r 180x360  remap to lat/lon uni grid')
        print (name,' -r 181x360  remap to lat/lon cap grid')
        print (name,' -s scriptfile')
        print (name,' -t timeindex starting at 1 [0=default-last frame. -1=all times]')
        print (name,' -y ngl,mpl')



        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i"):
            inputfile = arg
        elif opt in ("-j"):
            inputfile2 = arg
        elif opt in ("-t"):
            timeindex=int(arg)
            timeindex=timeindex-1  # convert to zero-indexing
            # 1-index  zero-index: -2 = last data
            #  -1        -2          all data 
            #  0         -1          last frame
            #  1..N      0..N-1      specific fame
        elif opt in ("-k"):
            levindex=int(arg)
            if levindex>0:
                levindex=levindex -1   # convert to zero-indexing
        elif opt in ("-p"):
            pressurelev=numpy.array([float(arg)])
        elif opt in ("-c"):
            clev=[float(i) for i in arg.split(",")]
        elif opt in ("-r"):
            nlatlon_interp=[int(i) for i in arg.split("x")]
        elif opt in ("-m"):
            projection=arg
        elif opt in ("-y"):
            if arg == "mpl": use_ngl=False
        elif opt in ("-s"):
            scripfile = arg
        elif opt in ("-g"):
            gllfile = arg
        elif opt in ("-e"):
            se_file = arg
        elif opt in ("-f"):
            contour_opt = arg
        elif opt in ("-o"):
            coutlines = int(arg)

    print("inputfile=",inputfile)
    return inputfile,inputfile2,args,projection,timeindex,levindex,pressurelev,clev,\
        nlatlon_interp,use_ngl,scripfile,gllfile,se_file,contour_opt,coutlines


def interp_to_latlon(data2d,lat,lon,lat_i,lon_i):
    # interpolating in lat/lon space has issues. interpolate in
    # stereographic projection:
    #
    # input:
    #    unstructured 1D data:  data(ncol),lat(ncol),lon(ncol)
    #    target lat/lon grid:   lon_i(nlon), lat_i(nlat)
    #
    # output 2D interpolated data::
    #   data(nlon,nlat)
    #
    
    dproj=crs.PlateCarree()
    halo = 15 # degrees        take hemisphere + halo for source grid

    if  lat_i[0]<0 and lat_i[-1]>0:
        # split grid into NH and SH
        # mesh grid
        nhalf = int(len(lat_i)/2)
        lat_south = lat_i[ :nhalf]
        lat_north = lat_i[ nhalf:]
        
        # take source data in the correct hemisphere, include extra halo points for interpolation
        # using the full global data sometimes confuses griddata with points being mapped close to infinity
        data2d_h=data2d[lat<halo]
        
        lon_h=lon[lat<halo]
        lat_h=lat[lat<halo]
        xv,yv=numpy.meshgrid(lon_i,lat_south)
        coords_in  = crs.SouthPolarStereo().transform_points(dproj,lon_h,lat_h)
        coords_out = crs.SouthPolarStereo().transform_points(dproj,xv.flatten(),yv.flatten())
        data_s = griddata(coords_in[:,0:2], data2d_h, coords_out[:,0:2], method='linear')
        
        data2d_h=data2d[lat>-halo]
        lon_h=lon[lat>-halo]
        lat_h=lat[lat>-halo]
        xv,yv=numpy.meshgrid(lon_i,lat_north)
        coords_in  = crs.NorthPolarStereo().transform_points(dproj,lon_h,lat_h)
        coords_out = crs.NorthPolarStereo().transform_points(dproj,xv.flatten(),yv.flatten())
        data_n = griddata(coords_in[:,0:2], data2d_h, coords_out[:,0:2], method='linear')
        
        data_i=numpy.concatenate((data_s,data_n)).reshape(len(lat_i),len(lon_i))

    elif lat_i[-1]<0:
        # SH only
        data2d_h=data2d[lat<halo]
        lon_h=lon[lat<halo]
        lat_h=lat[lat<halo]
        
        xv,yv=numpy.meshgrid(lon_i,lat_i)
        coords_in  = crs.SouthPolarStereo().transform_points(dproj,lon_h,lat_h)
        coords_out = crs.SouthPolarStereo().transform_points(dproj,xv.flatten(),yv.flatten())
        data_i = griddata(coords_in[:,0:2], data2d_h, coords_out[:,0:2], method='linear')

    elif lat_i[0]>0:
        # NH only
        data2d_h=data2d[lat>-halo]
        lon_h=lon[lat>-halo]
        lat_h=lat[lat>-halo]
        
        xv,yv=numpy.meshgrid(lon_i,lat_i)
        coords_in  = crs.NorthPolarStereo().transform_points(dproj,lon_h,lat_h)
        coords_out = crs.NorthPolarStereo().transform_points(dproj,xv.flatten(),yv.flatten())
        data_i = griddata(coords_in[:,0:2], data2d_h, coords_out[:,0:2], method='linear')

    else:
        print("Error: interp_to_latlon failed")
        sys.exit(1)
        
    return data_i
    

def extract_level(dataf,klev,plev,PS,hyam,hybm,kidx=0):

    if kidx==0:  # intepolate first index:
        if plev == None:
            data2d=dataf[klev,...]
        else:
            print("Interpolating 1st dimension")
            # vertical interpolation
            v_interp = 2   # type of interpolation: 1 = linear, 2 = log, 3 = loglog
            extrap = True  # is extrapolation desired if data is outside the range of PS
            P0mb = 1000    # ps in Pa, but this argument must be in mb
            
            if len(dataf.shape)==2:   # lev,ncol
                dataf2=numpy.expand_dims(dataf,axis=2)  # lev,ncol,1
                PS2=numpy.expand_dims(PS,axis=1)        # ncol,1
                data2d=numpy.squeeze(Ngl.vinth2p(dataf2,hyam,hybm,plev,PS2,v_interp,P0mb,1,extrap))
            elif len(dataf.shape)==3:  # lev,lat,lon
                data2d=numpy.squeeze(Ngl.vinth2p(dataf,hyam,hybm,plev,PS,v_interp,P0mb,1,extrap))
            else:
                print("ERROR: extract_level: dataf() needs to be 2 or 3 dimensiosn")
                
        return data2d

    if kidx==len(dataf.shape)-1:
        print("need to interpolating last dimension!")
        if plev == None:
            data2d=dataf[...,klev]
        else:
            # vertical interpolation
            print("vertical interpolation last dimension, not coded")
            sys.exit(1)
        return data2d

    print("Error: level dimension was not first or last")
    sys.exit(1)


def ngl_plot(wks,data2d,lon,lat,title,longname,units,
             projection,clev,cmap,scrip_file,se_file,contour_opt,coutlines,data2d_2=numpy.array([])):

    se_num=0
    if os.path.isfile(se_file):
        print("adding to plot:  SE grid from:",se_file)
        infile = Dataset(se_file,"r")
        se_coord  = infile.variables["coord"][:,:]
        se_connect  = infile.variables["connect1"][:,:]
        se_num=se_connect.shape[0]
        print("number of elements =",se_num)
        #print("se_coord shape: ", se_coord.shape)
        #print("se_connectshape: ", se_connect.shape)
        # j=se_connect(i,0:3) is the index of the 4 corners of cell i
        # se_coord(0:2,j) are the x,y,z coords of vertex j
        se_lat = numpy.arcsin(se_coord[2,:])*180/numpy.pi
        se_lon = numpy.arctan2(se_coord[1,:],se_coord[0,:])*180/numpy.pi
        

    cellbounds=False
    if os.path.isfile(scrip_file):
        infile = Dataset(scrip_file,"r")
        clat  = infile.variables["grid_corner_lat"][:,:]
        clon  = infile.variables["grid_corner_lon"][:,:]
        if clon.shape[0] == len(lat):
            cellbounds=True
    
    res = Ngl.Resources()

    #if cmap!=None:
    res.cnFillPalette   = cmap

    nlevels=-1  # not specified
    if len(clev)==1:
        res.cnLevelSelectionMode  = "AutomaticLevels"
        res.cnMaxLevelCount = clev[0]
        nlevels=clev[0]
    elif len(clev)==4:
        # log spacing
        cbase=clev[3]
        c0=numpy.log(clev[0])/numpy.log(cbase)
        c1=numpy.log(clev[1])/numpy.log(cbase)
        cinc=clev[2]

        #nlevels=(clev[1]-clev[0])/clev[2]
        nlevels=(c1-c0)/cinc
        clevs=[ c0+i*cinc  for i in range(1+round(nlevels))]
        clevs=numpy.power(cbase,clevs)
        print("nlevels=",nlevels," lev spacing ratio= ",cbase,"^",cinc)
        print("log clevs=",clevs)

        res.cnLevelSelectionMode = "ExplicitLevels" 
        res.cnLevels=clevs
    else:
        res.cnLevelSelectionMode  = "ManualLevels"
        res.cnMinLevelValF=clev[0]
        res.cnMaxLevelValF=clev[1]
        res.cnLevelSpacingF=clev[2]
        nlevels=(clev[1]-clev[0])/clev[2]




# defaults. some projection options might change:
    res.cnFillOn              = True           # Turn on contour fill.
    if contour_opt=='lo':
        res.cnFillOn              = False      # Turn on contour fill off.

    res.cnLinesOn             = False          # Turn off contour lines
    if contour_opt=='lo' or contour_opt=='la':
        res.cnLinesOn             = True       # Turn on contour lines

    res.cnLineLabelsOn        = False          # Turn off line labels.
    res.cnInfoLabelOn         = False          # Turn off info label.

    res.mpOutlineOn          = True
    if coutlines==0:     res.mpOutlineOn          = False

    res.mpFillOn             = False
    res.mpGridAndLimbOn      = False    # dont draw grid lines
    #res.mpShapeMode          = "FreeAspect"
    



    if projection == "latlon" :
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "MaximalArea"
        #    res.mpMinLatF = -90.
        #    res.mpMaxLatF = 90.
        #    res.mpMinLonF = -180.
        #    res.mpMaxLonF = 180.
    elif projection == "latlon-nc1":
        # no continenents, turn on contour lines
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "MaximalArea"
    elif projection == "latlon-nc2":
        # no continenents, turn off contour lines
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "MaximalArea"
    elif projection == "US1":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -30
        res.mpMaxLatF = 75
        res.mpMinLonF = -180
        res.mpMaxLonF =  0
    elif projection == "europe":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = 20.
        res.mpMaxLatF = 75.
        res.mpMinLonF = -40.
        res.mpMaxLonF =  40.
    elif projection == "andes":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -40.
        res.mpMaxLatF = 15.
        res.mpMinLonF = -100.
        res.mpMaxLonF =  -40.
    elif projection == "pacific1":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -20.
        res.mpMaxLatF = 60.
        res.mpMinLonF = -180.
        res.mpMaxLonF =  -80.
    elif projection == "andes2":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -17.
        res.mpMaxLatF =  3.
        res.mpMinLonF = -85.
        res.mpMaxLonF =  -65.
    elif projection == "himalaya":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = 0.
        res.mpMaxLatF = 60.
        res.mpMinLonF = 50.
        res.mpMaxLonF = 110.
    elif projection == "oro":
        res.mpProjection      = "Orthographic"   
        res.mpCenterLatF      =  45.
        res.mpCenterLonF      = -45.
    elif projection == "oro-debug1":
        res.mpProjection      = "Orthographic"   
        res.mpLimitMode="Angles"
        hdeg=10
        res.mpLeftAngleF=hdeg
        res.mpRightAngleF=hdeg
        res.mpBottomAngleF=hdeg
        res.mpTopAngleF=hdeg
    elif projection == "debug1":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = 30.
        res.mpMaxLatF = 50.
        res.mpMinLonF = 85.
        res.mpMaxLonF = 105.
    elif projection == "debug2":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -10.
        res.mpMaxLatF = 75.
        res.mpMinLonF = 45.
        res.mpMaxLonF = 175.
    elif projection == "debug3":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = 10.
        res.mpMaxLatF = 50.
        res.mpMinLonF = 60.
        res.mpMaxLonF = 100.
    elif projection == "baroclinic":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpCenterLonF         = 100.
        res.mpMinLatF = 25.
        res.mpMaxLatF = 75.
        res.mpMinLonF = 25.
        res.mpMaxLonF = 175.
    elif projection == "barotopo":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        #res.mpCenterLonF         = -90.
        res.mpMinLatF = 15. 
        res.mpMaxLatF = 75. 
        res.mpMinLonF = -150.
        res.mpMaxLonF = 20.
    elif projection == "barotopo-nc2":
        # no continenents, turn off contour lines
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        #res.mpCenterLonF         = -90.
        res.mpMinLatF = 15. 
        res.mpMaxLatF = 75. 
        res.mpMinLonF = -150.
        res.mpMaxLonF = 20.
    else:
        print("Bad projection argument: ",projection)
        sys.exit(3)

        
    res.nglFrame = False # Don't advance frame.
    #res.tiXAxisString = "~F25~longitude"
    #res.tiYAxisString = "~F25~latitude"
    res.nglPointTickmarksOutward = True
    

    if cellbounds:
        res.cnFillMode = 'CellFill'
        res.sfXCellBounds = clon
        res.sfYCellBounds = clat
    else:
        if contour_opt=='area' or contour_opt=='la':
            res.cnFillMode            = "AreaFill"
        elif contour_opt=='raster':
            res.cnFillMode            = "RasterFill"
        else:
            res.cnFillMode            = "AreaFill"
            if (nlevels>25):
                print("more than 25 contours, switching to raserfill")
                res.cnFillMode            = "RasterFill"
        res.cnRasterSmoothingOn = True
        res.sfXArray = lon[:]
        res.sfYArray = lat[:]



    #res.sfCopyData = False
 
#    "not a valid resource in contour at this time...       
#    if wks_type == "pdf":
#        res.gsnMaximize           = True        
#        res.gsnPaperOrientation   = "portrait"

    #res.lbLabelAutoStride   = True         # Clean up labelbar labels.
    #res.lbAutoManage = True
    #res.lbLabelStride       = 10
    res.lbBoxLinesOn        = False        # Turn of labelbar box lines.
    res.lbOrientation       = "horizontal"
    


    print("Title: ",title)
    print("Longname: ",longname)
    
    res.tiMainString = title

    print("data min/max=",numpy.amin(data2d),numpy.amax(data2d))        
    if res.cnLevelSelectionMode == "ManualLevels":
        print("contour levels: manual [",res.cnMinLevelValF,",",\
              res.cnMaxLevelValF,"] spacing=",res.cnLevelSpacingF)
        print("number of contour levels:",nlevels)
    elif res.cnLevelSelectionMode=="AutomaticLevels":
        print("contour levels: auto. number of levels:",res.cnMaxLevelCount)


    if nlevels>20:
        res.lbLabelStride       = nlevels/8
        

    # for lat/lon plots, add cyclic point:
    res2=res
    if hasattr(res,"sfXArray"):
        if len(res.sfXArray)*len(res.sfYArray) == numpy.prod(data2d.shape):
            print("NGL structured data plot.  Adding cyclic point")
            data2d,lon2= Ngl.add_cyclic(data2d[:,:],res.sfXArray[:])
            res2.sfXArray=lon2
        else:
            print("NGL unstructered plot with internal triangulation")
    elif hasattr(res,"sfXCellBounds"):
        print("NGL unstructered plot with cell bounds")
    else:
        print("Error with resource coordinate data")

    if data2d_2.size == 0:
        # contour data2d
        map1 = Ngl.contour_map(wks,data2d,res2)
    else:
        # plot variable, then add contour lines from data2d_2:
        res2.nglDraw  = False
        map1 = Ngl.map(wks,res2)
        map2 = Ngl.contour(wks,data2d,res2)
        print("Adding contour line plot from data2d_2...")
        print("data2d_2 min/max=",numpy.amin(data2d_2),numpy.amax(data2d_2))

        res3=res
        res3.mpOutlineOn          = False
        res3.cnFillOn             = False         # Turn off contour fill.
        res3.cnLinesOn            = True          # Turn on contour lines
        #res3.cnLineColor          = "White"

        #res3.cnLevelSelectionMode  = "AutomaticLevels"
        #res3.cnMaxLevelCount = 5
        #res3.cnLevelSpacingF=1e5   # ignored, but set to prevent warnings

        res3.cnMinLevelValF=200
        res3.cnMaxLevelValF=300
        res3.cnMaxLevelCount = 100
        res3.cnLevelSpacingF=5  
        res3.cnLineThicknessF = 1.5

        map3 = Ngl.contour(wks,data2d_2,res3)
        Ngl.overlay(map1,map2)
        Ngl.overlay(map1,map3)
        Ngl.maximize_plot(wks,map1)    # Maximize size of plot in frame.
        Ngl.draw(map1)
        del res3

    del res2
    print("Contour done.")
        
    #-- write variable long_name and units to the plot
    txres = Ngl.Resources()
    txres.txFontHeightF = 0.016
    Ngl.text_ndc(wks,longname,0.20,0.95,txres)
    Ngl.text_ndc(wks,units, 0.85,0.95,txres)
    del txres

    if se_num>0:
        gsres = Ngl.Resources()
        gsres.gsLineColor            = "black"
        gsres.gsLineThicknessF       = 0.05
        for x in se_connect:
            x1=se_lon[x-1]
            x2=x1[0:1]
            plon=numpy.concatenate( [x1,x2])
            x1=se_lat[x-1]
            x2=x1[0:1]
            plat=numpy.concatenate( [x1,x2])

            # se_coord[0:2,x-1] = for element corners in cartesian coordinate.
            # compute diagonal distance:
            #xcoord=se_coord[:,x-1]
            #d1=numpy.sum((xcoord[:,0]-xcoord[:,2])**2)
            #d2=numpy.sum((xcoord[:,1]-xcoord[:,3])**2)
            #dist2 = ((numpy.sqrt(d1/2)+numpy.sqrt(d2/2))/2)* 6.4e3  # dist in km
            #if (dist2> 80):    # NE30 elements are about 325km
            Ngl.polyline(wks,map1,plon,plat,gsres)
        del gsres
            
    
    Ngl.frame(wks)       # advance frame
    return map1


    

def mpl_plot(data2d,lon,lat,title,longname,units,proj,clev,cmap,scrip_file,gllfile,contour_opt,coutlines):
    
    # Setup the plot
    figure = pyplot.figure() #(figsize=(15, 10))
    dataproj=crs.PlateCarree()

    # pcolor/tripcolor doesn't use nelvels or contour intervals
    if len(clev)==1:
        vmin=None
        vmax=None
        nlevels=int(round(clev[0]))
        levels=nlevels  # matplotlib will choose values
    elif len(clev)==3:
        vmin=clev[0]
        vmax=clev[1]
        nlevels=int(round( (clev[1]-clev[0])/clev[2] ))
        levels=[vmin+i*clev[2] for i in range(nlevels+1)]
        print("levels=",levels)
    #elif len(clev)==4:
    #not all MPL plotting options support list of levels
    else:
        print("Error: mpl len(clev) <> 1,3")
        sys.exit(1)

            
    print("num contour levels=",nlevels)

    # set_extent[lon_min,lon_max,lat_min,lat_mx]
    if proj=="latlon":
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_global()
    elif proj=="US1":
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([-180, 0, -30, 75],crs=dataproj)
    elif proj=="andes":
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([-100, -40, -40, 15],crs=dataproj)
    elif proj=="himalaya":
        plotproj=crs.PlateCarree(central_longitude=90.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([50, 110, 0, 60],crs=dataproj)
    elif proj=="oro":
        plotproj=crs.Orthographic(central_longitude=-45.0, central_latitude=45.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_global()
    elif proj == "europe":
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([-40, 40, 20, 75],crs=dataproj)
    elif proj == "debug3":
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([60, 100, 10, 50],crs=dataproj)
    else:
        print("Bad projection argument: ",projection)
        sys.exit(3)

    # tried to download data set        
    #ax.coastlines(linewidth=0.2)
    
    # structured lat/lon or unstructured data?
    struct=False
    if len(lon)*len(lat) == numpy.prod(data2d.shape): struct=True
    
    compute_tri=True
    if ~struct and os.path.isfile(gllfile):
        cfile = Dataset(gllfile,"r")
        ec=cfile.variables["element_corners"]
        nd=ec.shape
        # by Euler, number of subcells is number of gll nodes -2
        if (nd[1] == len(lat)-2):
            ntris=nd[1]*2
            tri=numpy.empty((ntris,3), dtype=int)
            
            tri[::2,0]=ec[0,:]
            tri[::2,1]=ec[1,:]
            tri[::2,2]=ec[2,:]
            
            tri[1::2,0]=ec[0,:]
            tri[1::2,1]=ec[2,:]
            tri[1::2,2]=ec[3,:]
            tri=tri-1   # zero indexing
            compute_tri=False

    cellbounds=False
    if os.path.isfile(scrip_file):
        infile = Dataset(scrip_file,"r")
        clat  = infile.variables["grid_corner_lat"][:,:]
        clon  = infile.variables["grid_corner_lon"][:,:]
        # if in Radians, convert to degrees
        if "radian" in infile.variables["grid_corner_lat"].units.lower():
            clon=clon*180/numpy.pi
            clat=clat*180/numpy.pi
        if clon.shape[0] == len(lat):
            cellbounds=True
            compute_tri=False

    pl2=None
    print("colormap min/max=",vmin,vmax)
    print("data min/max=",numpy.amin(data2d),numpy.amax(data2d))
    if struct:
        data2d_ext, lon2 = add_cyclic_point(data2d, coord=lon,axis=1)
        print("MPL plotting structured data (with added cyclic point)")
        if contour_opt=='' or contour_opt=='raster':
            print("using pcolormesh")
            pl=ax.pcolormesh(lon2, lat, data2d_ext,levels,vmin=vmin,vmax=vmax,
                             transform=dataproj, cmap=cmap)
        elif contour_opt=='la':
            print("contour lines + area fill")
            pl=ax.contourf(lon2, lat, data2d_ext, levels,vmin=vmin,vmax=vmax,
                           transform=dataproj, cmap=cmap)
            pl2=ax.contour(pl, levels,vmin=vmin,vmax=vmax,transform=dataproj,
                           colors='k',linewidths=.5)
        elif contour_opt=='lo':
            print("contour lines only")
            pl=ax.contour(lon2, lat, data2d_ext, levels,vmin=vmin,vmax=vmax,
                           transform=dataproj, colors='k',linewidths=.5)
        elif contour_opt=='area':
            print("using contourf (fill only)")
            pl=ax.contourf(lon2, lat, data2d_ext, levels,vmin=vmin,vmax=vmax,
                           transform=dataproj, cmap=cmap)
        
    elif compute_tri:
        print("MPL plot using internal Delaunay triangulation")
        # do the triangulation in the plot coordinates for better results
        tcoords = plotproj.transform_points(dataproj,lon[:],lat[:])
        # need to remove non-visible points
        xi=tcoords[:,0]!=numpy.inf
        tc=tcoords[xi,:]
        datai=data2d[:][xi]  # convert to numpy array, then subset
        
        pl = ax.tripcolor(tc[:,0],tc[:,1], datai,vmin=vmin, vmax=vmax,
                          shading='gouraud',cmap=cmap)

        #
        # as with above, we could put in options for tricontour & tricontourf
        #
    elif cellbounds:
        from matplotlib.collections import PolyCollection
        print("MPL plot using scrip cells")
        # latlon->cartesian->local coords. this will put any seams at plot boundaries
        proj3d=crs.Geocentric()   # for cartesian (x,y,z) representation
        x3d = proj3d.transform_points(dataproj,clon[:,:],clat[:,:])
        # x3d[:,:,0:3] x,y,z coords of 
        ccoords = plotproj.transform_points(proj3d,x3d[:,:,0],x3d[:,:,1],x3d[:,:,2])
        # ccoords will be [:,:,0:2]  for lat,lon,r but r=0
        # remove bad cells:
        # cells that stradle any coordinate branch cut, or
        # are non-visible
        x0=numpy.amin(ccoords[:,:,:],1)
        x1=numpy.amax(ccoords[:,:,:],1)
        d= ((x0[:,0]-x1[:,0])**2 + (x0[:,1]-x1[:,1])**2)**0.5
        dmax=max(d)
        dmin=min(d)
        print("cell diameter min,max:",dmin,dmax)
        xi=(d < 100*dmin)  # and (x1!=numpy.inf) and (y1!=numpy.inf)
        datai=data2d[:][xi]
        print("cells removed: ",len(data2d)-len(datai),"out of",len(data2d))
        corners=numpy.stack([ccoords[xi,:,0],ccoords[xi,:,1]],axis=2)
        
        # create cells
        # antialized=False needed to avoid visable cell edges
        # but when combined with alpha/transparancy, cell edges become visiable again
        # due to some issue with how alpha is applied (or not applied) at edges
        p = PolyCollection(corners,array=datai,edgecolor='none',linewidths=0,antialiased=False)
        p.set_clim([vmin,vmax])
        pl=ax.add_collection(p)
    else:

        print("MPL plot using gll subcell triangulation")
        # latlon->cartesian->local coords. this will put any seams at plot boundaries
        proj3d=crs.Geocentric()   # for cartesian (x,y,z) representation
        x3d = proj3d.transform_points(dataproj,lon[:],lat[:])
        tcoords = plotproj.transform_points(proj3d,x3d[:,0],x3d[:,1],x3d[:,2])

        #Remove bad triangles:
        x0=tcoords[tri[:,0],0]
        y0=tcoords[tri[:,0],1]
        x1=tcoords[tri[:,1],0]
        y1=tcoords[tri[:,1],1]
        x2=tcoords[tri[:,2],0]
        y2=tcoords[tri[:,2],1]
        d=numpy.empty(tri.shape)
        d[:,0]=((x0-x1)**2 + (y0-y1)**2)**0.5
        d[:,1]=((x0-x2)**2 + (y0-y2)**2)**0.5
        d[:,2]=((x1-x2)**2 + (y1-y2)**2)**0.5
        dmax=numpy.amax(d,axis=1)
        gmin=numpy.nanmin(dmax[dmax != numpy.inf])
        gmax=numpy.nanmax(dmax[dmax != numpy.inf])
        print("triangle max lengths: ",gmin,gmax)
        mask = numpy.logical_or( dmax > 25*gmin, numpy.isnan(dmax))
        # gouraud shading requires we remove non-visable triangles
        pl = ax.tripcolor(tcoords[:,0],tcoords[:,1],tri,data2d,vmin=vmin,vmax=vmax,
                          mask=mask,shading='gouraud',cmap=cmap)
        # plot some of the triangles to make sure they are ok:
        #ax.triplot(tcoords[:,0],tcoords[:,1],tri[1:100,:],'go-')

        # if we specify the trianglization "tri", then can also give data per vertex or per face
        # for per-face data, might need facecolors=zfaces option.  but for facedata, probably 
        # faster to just use the above polycollection, which should be identical
        # HOWEVER: tripcolor+facecolors may avoid the edge effects with Polycollection + transparancey

    pl.set_clim([vmin,vmax])
    if units=="":
        label='%s'%(longname)
    else:
        label='%s (%s)'%(longname, units)

    cb = pyplot.colorbar(pl, orientation='horizontal', 
                        label=label,shrink=0.75, pad=0.1)
    # add in contour lines to the color bar:
    #if pl2!=None: cb.add_lines(pl2)

    # Plot GLL nodes for perspective
    #pl = ax.plot(lon, lat, 'k.', projection=dataproj, markersize=1)



def mpl_streamlines(data2d,data2d_2,lon,lat,title,longname,units,proj,clev,cmap):
    
    # Setup the plot
    figure = pyplot.figure() #(figsize=(15, 10))
    dataproj=crs.PlateCarree()

    # pcolor/tripcolor doesn't use nelvels or contour intervals
    if len(clev)==1:
        vmin=None
        vmax=None
        nlevels=int(round(clev[0]))
    else:
        vmin=clev[0]
        vmax=clev[1]
        nlevels=int(round( (clev[1]-clev[0])/clev[2] ))
    print("num contour levels=",nlevels)

    # set_extent[lon_min,lon_max,lat_min,lat_mx]
    if proj=="latlon":
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_global()
    elif proj=="US1":
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([-180, 0, -30, 75],crs=dataproj)
    elif proj=="andes":
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([-100, -40, -40, 15],crs=dataproj)
    elif proj=="himalaya":
        plotproj=crs.PlateCarree(central_longitude=90.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([50, 110, 0, 60],crs=dataproj)
    elif proj=="oro":
        plotproj=crs.Orthographic(central_longitude=-45.0, central_latitude=45.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_global()
    elif proj == "europe":
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([-40, 40, 20, 75],crs=dataproj)
    elif proj == "debug3":
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([60, 100, 10, 50],crs=dataproj)
    else:
        print("Bad projection argument: ",proj)
        sys.exit(3)

    # structured lat/lon or unstructured data?
    struct=False
    if len(lon)*len(lat) == numpy.prod(data2d.shape): struct=True
    if not struct:
        print("streamlines requires lat/lon data")
        sys.exit(2)
        
    print("colormap min/max=",vmin,vmax)
    print("data min/max=",numpy.amin(data2d),numpy.amax(data2d))

    
    data2d_ext, lon2 = add_cyclic_point(data2d, coord=lon,axis=1)
    data2d_2_ext, lon3 = add_cyclic_point(data2d_2, coord=lon,axis=1)
    magnitude = numpy.sqrt(data2d_ext**2 + data2d_2_ext**2)
    pl=ax.streamplot(lon2, lat, data2d_ext, data2d_2_ext, transform=dataproj,  linewidth=2, density=2, color=magnitude)

    if units=="":
        cb = pyplot.colorbar(pl.lines, orientation='horizontal', 
                             label='%s'%(longname),shrink=0.75, pad=0.1)
    else:
        cb = pyplot.colorbar(pl.lines, orientation='horizontal', 
                             label='%s (%s)'%(longname, units),shrink=0.75, pad=0.1)


