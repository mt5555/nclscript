import numpy, os, sys, getopt

# needed for mpl_plot
from cartopy import crs
from cartopy.util import add_cyclic_point
import matplotlib.tri as mpltri
from matplotlib import pyplot
from scipy.interpolate import griddata
# needed for ngl_plot
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
    coutlines=0
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
    print("colormap min/max=",vmin,vmax)

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
    if struct:
        data2d_ext, lon2 = add_cyclic_point(data2d, coord=lon,axis=1)
        print("MPL plotting structured data (with added cyclic point)")
        if contour_opt=='' or contour_opt=='raster':
            print("using pcolormesh")
            pl=ax.pcolormesh(lon2, lat, data2d_ext,vmin=vmin,vmax=vmax,
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

        if ("proj=eqc" in plotproj.srs) or ("proj=robin" in plotproj.srs):
            # all polygons always visable
            tc=tcoords
            datai=data2d
        else:   # if "proj=ortho" in plotproj.srs:
            # Triangulation chokes if there are  non-visible points
            xi=numpy.logical_and ( numpy.isfinite(tcoords[:,0]), numpy.isfinite(tcoords[:,1]))
            tc=tcoords[xi,:]
            datai=data2d[xi]  


        # compute triangularization
        tri=mpltri.Triangulation(tc[:,0],tc[:,1])

        # diam = numpy.hypot(tc[tri.triangles,0].mean(axis=1),
        #                 tc[tri.triangles,1].mean(axis=1))
        #gmin=numpy.nanmin(diam[diam != numpy.inf])
        #gmax=numpy.nanmax(diam[diam != numpy.inf])
        #print("triangle max lengths: ",gmin,gmax)
        #tri.set_mask( numpy.logical_or(diam > 95*gmin, numpy.isnan(diam)))

        #tripcolor requires vmin/vmax set

        # gouraud shading requires we remove non-visable triangles        
        # gouraud shading seems broken for pdf, be sure to use png
        if contour_opt=='' or contour_opt=='raster':
            if vmin==None:
                vmin=numpy.amin(data2d)
                vmax=numpy.amax(data2d)
                print("using data min/max for colormap, min/max=",vmin,vmax)
            print("using tripcolor")
            pl = ax.tripcolor(tri, datai,vmin=vmin, vmax=vmax,
                              shading='gouraud',cmap=cmap,antialiased=False)
        elif contour_opt=='la':
            print("contour lines + area fill")
            pl = ax.tricontourf(tri, datai,levels, vmin=vmin, vmax=vmax,cmap=cmap)
            pl = ax.tricontour(tri, datai,levels, vmin=vmin, vmax=vmax, 
                                colors='k',linewidths=.5)
        elif contour_opt=='lo':
            print("contour lines only")
            pl = ax.tricontour(tri, datai,levels, vmin=vmin, vmax=vmax,
                                colors='k',linewidths=.5)
        elif contour_opt=='area':
            print("using contourf (fill only)")
            pl = ax.tricontourf(tri, datai,levels, vmin=vmin, vmax=vmax,cmap=cmap)

        #ax.triplot(tri,'g-',lw=.1) # plot triangles for debugging
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

        # better approach: see ../polycollection/plotpoly_mpl.py
        # for lat/lon & robin, look for cut cells, fix them to +180, and duplicate at -180
        # for other projections, just mask out points with numpy.isfinite=False
        # no need to project to 3D first
        #


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
        datai=data2d[xi]
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

        # better approach: see ../polycollection/plotpoly_mpl.py
        # for lat/lon & robin, look for cut cells, fix them to +180, and duplicate at -180
        # for other projections, just mask out points with numpy.isfinite=False
        # no need to project to 3D first
        #

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


