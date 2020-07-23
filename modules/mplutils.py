import numpy, os, sys, getopt
from cartopy import crs
from cartopy.util import add_cyclic_point
from matplotlib import pyplot

#
# setup resources for contour plots
#
# can handle 3 types of plots:
#
# 1. latlon data:  data(nlat,nlon) with lat(nlat) and lon(nlon)
# 2. unstructured:  data(ncol), lat(ncol), lon(ncol)
# 3. unstructured with cell fill:  data(ncol), lat(ncol), lon(ncol)
#        and scrip_file is a file containing the cell bounds:
#        grid_corner_lat(ncol,4), grid_corner_lon(ncol,4)
#
#

def myargs(argv):
    inputfile = ''
    scripfile = ''
    gllfile = ''
    varname = ''
    name = argv[0]
    try:
        opts, args = getopt.getopt(argv[1:],"i:s:g:")
    except getopt.GetoptError:
        print (name,' -i inputfile [-s scriptfile] [-g gll_subcell_file ]  varname')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i"):
            inputfile = arg
        elif opt in ("-s"):
            scripfile = arg
        elif opt in ("-g"):
            gllfile = arg
                
    return inputfile,args,scripfile,gllfile


def mpl_plot(data2d,lon,lat,proj,nlevels,cmap,gllfile):
    
    # Setup the plot
    figure = pyplot.figure(figsize=(15, 10))
    dataproj=crs.PlateCarree()
    

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
        print("Plotting structured data")
        data2d_ext, lon2 = add_cyclic_point(data2d[:,:], coord=lon[:],axis=1)
        pl=ax.contourf(lon2, lat, data2d_ext, nlevels,transform=dataproj, cmap=cmap)
    elif compute_tri:
        print("Plot using internal Delaunay triangulation")
        # do the triangulation in the plot coordinates for better results
        tcoords = plotproj.transform_points(dataproj,lon[:],lat[:])
        # need to remove non-visible points
        xi=tcoords[:,0]!=numpy.inf
        tc=tcoords[xi,:]
        datai=data2d[:][xi]  # convert to numpy array, then subset
        
        pl = ax.tripcolor(tc[:,0],tc[:,1], datai,shading='gouraud',cmap=cmap)
    else:
        print("Using gll subcell triangulation")
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
        pl = ax.tripcolor(tcoords[:,0],tcoords[:,1],tri,data2d,mask=mask,shading='gouraud',cmap=cmap)
        # plot some of the triangles to make sure they are ok:
        #ax.triplot(tcoords[:,0],tcoords[:,1],tri[1:100,:],'go-')
        
    
    longname=""
    units=""
    if hasattr(data2d,"longname"):
        longname=data2d.longname
    if hasattr(data2d,"units"):
        units=data2d.units
    cb = pyplot.colorbar(pl, orientation='horizontal', 
                         label='%s (%s)'%(longname, units),shrink=0.75, pad=0.1)

    # Plot GLL nodes for perspective
    #pl = ax.plot(lon, lat, 'k.', projection=dataproj, markersize=1)
    
