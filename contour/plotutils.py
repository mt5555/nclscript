import numpy, os, sys, getopt

# needed for mpl_plot
from cartopy import crs
from cartopy.util import add_cyclic_point
from matplotlib import pyplot
from scipy.interpolate import griddata

# needed for ngl_plot
import Ngl
import Nio

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
# 
#
# can handle many types of plots:
#
# 1. latlon data:  data(nlat,nlon)  NGL and MPL.   best results.
# 2. unstructured  dual grid:  NGL only
#                  2a: construct dual grid via triangulation.  UGLY
#                  2b: specify dual grid from SCRIP file.  BEST FOR PG2 data
# 3. unstructured: vertex data:  MPL only
#                  3a: construct triangulation in plot coordinates BEST FOR GLL DATA
#                  3b: construct triangulation from subcell grid
#
#

def myargs(argv):
    inputfile = ''
    scripfile = ''
    gllfile = ''
    varname = ''
    use_ngl = True
    timeindex = None
    levindex = None
    pressurelev = None
    nlatlon_interp=None
    clev = None
    name = argv[0]
    projection='latlon'
    try:
        opts, args = getopt.getopt(argv[1:],"i:s:g:t:k:p:y:c:m:r:")
    except getopt.GetoptError:
        print (name,' -i inputfile [options] varname')
        print (name,' -t timeindex starting at 1 [0=default-all times. -1=last time]')
        print (name,' -k levindex starting at 1  [default: 3*nlev/4]')
        print (name,' -p pressure(mb)  interpolate to pressure level')
        print (name,' -c nlevels  number of contour levels (ignored in MPL)')
        print (name,' -c cmin,cmax       contour level min,max with 40 levels')
        print (name,' -c cmin,cmax,cinc  contour level min,max,spacing')
        print (name,' -m map projeciton  latlon,US1,oro,andes,hamalaya,etc...')
        print (name,' -r 180x360  remap to lat/lon uni grid')
        print (name,' -r 181x360  remap to lat/lon cap grid')
        print (name,' -y ngl,mpl')
        print (name,' -s scriptfile')
        print (name,' -g gll_subcell_file')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i"):
            inputfile = arg
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
                
    return inputfile,args,projection,timeindex,levindex,pressurelev,clev,\
        nlatlon_interp,use_ngl,scripfile,gllfile


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
    
    # mesh grid
    dproj=crs.PlateCarree()
    nhalf = int(len(lat_i)/2)
    lat_south = lat_i[ :nhalf]
    lat_north = lat_i[ nhalf:]

    # take source data in the correct hemisphere, include extra halo points for interpolation
    # using the full global data sometimes confuses griddata with points being mapped close to infinity
    halo = 15 # degrees
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
    return data_i
    

def extract_level(dataf,klev,plev,PS,hyam,hybm):
    if plev == None:
        data2d=dataf[klev,...]
    else:
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


def ngl_plot(wks,data2d,lon,lat,title,longname,units,
             projection,clev,cmap,scrip_file):
    
    cellbounds=False
    if os.path.isfile(scrip_file):
        infile = Nio.open_file(scrip_file,"r")
        clat  = infile.variables["grid_corner_lat"][:,:]
        clon  = infile.variables["grid_corner_lon"][:,:]
        if clon.shape[0] == len(lat):
            cellbounds=True
    
    res = Ngl.Resources()

    if cmap!=None:
        res.cnFillPalette   = cmap
        
    if len(clev)==1:
        res.cnLevelSelectionMode  = "AutomaticLevels"
        res.cnMaxLevelCount = clev[0]
    else:
        res.cnLevelSelectionMode  = "ManualLevels"
        res.cnMinLevelValF=clev[0]
        res.cnMaxLevelValF=clev[1]
        res.cnLevelSpacingF=clev[2]

    if projection == "latlon" :
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "MaximalArea"
        #    res.mpMinLatF = -90.
        #    res.mpMaxLatF = 90.
        #    res.mpMinLonF = -180.
        #    res.mpMaxLonF = 180.
    elif projection == "US1":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -30
        res.mpMaxLatF = 75
        res.mpMinLonF = -180
        res.mpMaxLonF =  0
    elif projection == "andes":
        res.mpProjection = "CylindricalEquidistant"
        res.mpLimitMode = "LatLon"
        res.mpMinLatF = -40.
        res.mpMaxLatF = 15.
        res.mpMinLonF = -100.
        res.mpMaxLonF =  -40.
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
    else:
        print("Bad projection argument: ",projection)
        sys.exit(3)

        
    res.nglFrame = False # Don't advance frame.
    #res.tiXAxisString = "~F25~longitude"
    #res.tiYAxisString = "~F25~latitude"
    res.nglPointTickmarksOutward = True
    
    res.cnFillOn              = True           # Turn on contour fill.
    res.cnLinesOn             = False          # Turn off contour lines
    res.cnLineLabelsOn        = False          # Turn off line labels.
    res.cnInfoLabelOn         = False          # Turn off info label.

    if cellbounds:
        res.cnFillMode = 'CellFill'
        res.sfXCellBounds = clon
        res.sfYCellBounds = clat
    else:
        #res.cnFillMode            = "AreaFill"
        res.cnFillMode            = "RasterFill"
        res.cnRasterSmoothingOn = True
        res.sfXArray = lon[:]
        res.sfYArray = lat[:]

    #res.sfCopyData = False
 
#    "not a valid resource in contour at this time...       
#    if wks_type == "pdf":
#        res.gsnMaximize           = True        
#        res.gsnPaperOrientation   = "portrait"

    res.lbLabelAutoStride   = True         # Clean up labelbar labels.
    #res.lbAutoManage = True
    #res.lbLabelStride       = 10
    res.lbBoxLinesOn        = False        # Turn of labelbar box lines.
    res.lbOrientation       = "horizontal"
    
    res.mpOutlineOn          = True
    res.mpFillOn             = False
    res.mpGridAndLimbOn      = False    # dont draw grid lines
    #res.mpShapeMode          = "FreeAspect"
    


    
    res.tiMainString = title
    print("Title: ",res.tiMainString)

    print("data min/max=",numpy.amin(data2d),numpy.amax(data2d))        
    if res.cnLevelSelectionMode == "ManualLevels":
        nlevels=(res.cnMaxLevelValF-res.cnMinLevelValF)/res.cnLevelSpacingF
        print("contour levels: manual [",res.cnMinLevelValF,",",\
              res.cnMaxLevelValF,"] spacing=",res.cnLevelSpacingF)
        print("number of contour levels:",nlevels)
    else:
        print("contour levels: auto. number of levels:",res.cnMaxLevelCount)
        nlevels=res.cnMaxLevelCount

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
            
    # plot:
    map = Ngl.contour_map(wks,data2d,res2)
    print("Contour done.")
    del res2
        
    #-- write variable long_name and units to the plot
    txres = Ngl.Resources()
    txres.txFontHeightF = 0.016
    Ngl.text_ndc(wks,longname,0.14,0.82,txres)
    Ngl.text_ndc(wks,units, 0.95,0.82,txres)
    del txres
    
    Ngl.frame(wks)       # advance frame
    return map


    

def mpl_plot(data2d,lon,lat,title,longname,units,proj,clev,cmap,gllfile):
    
    # Setup the plot
    figure = pyplot.figure(figsize=(15, 10))
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

    if proj=="latlon":
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_global()
    elif proj=="US1":
        plotproj=plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([-180, 0, -30, 75],crs=dataproj)
    elif proj=="andes":
        plotproj=plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([-100, -40, -40, 15],crs=dataproj)
    elif proj=="himalaya":
        plotproj=plotproj=crs.PlateCarree(central_longitude=90.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_extent([50, 110, 0, 60],crs=dataproj)
    elif proj=="oro":
        plotproj=crs.Orthographic(central_longitude=-45.0, central_latitude=45.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_global()
    else:
        print("Bad projection argument: ",projection)
        sys.exit(3)

    ax.coastlines(linewidth=0.2)
    
    # strucgtured lat/lon or unstructured data?
    struct=False
    if len(lon)*len(lat) == numpy.prod(data2d.shape): struct=True
    
    compute_tri=True
    if ~struct and os.path.isfile(gllfile):
        cfile = Nio.open_file(gllfile,"r")
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


    print("data min/max=",numpy.amin(data2d),numpy.amax(data2d))
    print("colormap min/max=",vmin,vmax)
    if struct:
        data2d_ext, lon2 = add_cyclic_point(data2d, coord=lon,axis=1)
        print("MPL plotting structured data (with added cyclic point)")
        pl=ax.pcolormesh(lon2, lat, data2d_ext,vmin=vmin,vmax=vmax,
                       transform=dataproj, cmap=cmap)
        #pl=ax.contourf(lon2, lat, data2d_ext, nlevels,vmin=vmin,vmax=vmax,
        #               transform=dataproj, cmap=cmap)
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
        
    
    cb = pyplot.colorbar(pl, orientation='horizontal', 
                         label='%s (%s)'%(longname, units),shrink=0.75, pad=0.1)

    # Plot GLL nodes for perspective
    #pl = ax.plot(lon, lat, 'k.', projection=dataproj, markersize=1)
    
