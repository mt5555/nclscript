import numpy, os, sys, getopt

# needed for mpl_plot
from cartopy import crs
from cartopy.util import add_cyclic_point
from matplotlib import pyplot

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
    name = argv[0]
    try:
        opts, args = getopt.getopt(argv[1:],"i:s:g:t:")
    except getopt.GetoptError:
        print (name,' -i inputfile [-s scriptfile] [-g gll_subcell_file ] [-t ngl,mpl] varname')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i"):
            inputfile = arg
        elif opt in ("-t"):
            if arg == "mpl": use_ngl=False
        elif opt in ("-s"):
            scripfile = arg
        elif opt in ("-g"):
            gllfile = arg
                
    return inputfile,args,use_ngl,scripfile,gllfile



def ngl_plot(wks,wks_type,data2d,lon,lat,title,longname,units,
             projection,clev,cmap,scrip_file):
    
    cellbounds=False
    if os.path.isfile(scrip_file):
        infile = Nio.open_file(scrip_file,"r")
        clat  = infile.variables["grid_corner_lat"][:,:]
        clon  = infile.variables["grid_corner_lon"][:,:]
        if clon.shape[0] == len(lat):
            cellbounds=True
    
    res = Ngl.Resources()

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

        
    res.nglFrame = False # Don't advance frame.
    #res.tiXAxisString = "~F25~longitude"
    #res.tiYAxisString = "~F25~latitude"
    res.nglPointTickmarksOutward = True
    
    res.cnFillOn              = True           # Turn on contour fill.
    #res.cnFillPalette         = "BlueYellowRed"  # good for symmetric data
    #res.cnFillPalette         = 'MPL_viridis'
    #res.cnFillPalette         = "wgne15"
    #res.cnFillPalette         = "StepSeq25"
    res.cnFillPalette         = "WhiteBlueGreenYellowRed"
    #res.cnFillPalette         = "BlAqGrYeOrReVi200"
    res.cnLinesOn             = False          # Turn off contour lines
    res.cnLineLabelsOn        = False          # Turn off line labels.
    res.cnInfoLabelOn         = False          # Turn off info label.

    if cellbounds:
        res.cnFillMode = 'CellFill'
        res.sfXCellBounds = clon
        res.sfYCellBounds = clat
    else:
        res.cnFillMode            = "RasterFill"
        res.sfXArray = lon[:]
        res.sfYArray = lat[:]

    #res.sfCopyData = False
        
    if wks_type == "pdf":
        res.gsnMaximize           = True        
        res.gsnPaperOrientation   = "portrait"

    res.lbLabelAutoStride   = True         # Clean up labelbar labels.
    res.lbLabelStride       = 10 
    res.lbBoxLinesOn        = False        # Turn of labelbar box lines.
    res.lbOrientation       = "horizontal"
    
    res.mpOutlineOn          = True
    res.mpFillOn             = False
    res.mpGridAndLimbOn      = False    # dont draw grid lines
    #res.mpShapeMode          = "FreeAspect"
    


    
    res.tiMainString = title
    print("Title: ",res.tiMainString)


    if res.cnLevelSelectionMode == "ManualLevels":
        print("contour levels: manual [",res.cnMinLevelValF,",",\
              res.cnMaxLevelValF,"] spacing=",res.cnLevelSpacingF)
    else:
        print("contour levels: auto. number of levels:",res.cnMaxLevelCount)
        
    print("data min/max=",numpy.amin(data2d),numpy.amax(data2d))        

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

    if len(clev)==1:
        nlevels=clev[0].astype(int)
    else:
        vmin=clev[0]
        vmax=clev[1]
        nlevels=((clev[1]-clev[0])/clev[2] ).astype(int)

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
        print("Bad projection argument. assuming global lat/lon")
        plotproj=crs.PlateCarree(central_longitude=0.0)
        ax = pyplot.axes(projection=plotproj)
        ax.set_global()

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
    
