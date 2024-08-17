#
# hvscrip.py   test polygon plotting approaches from holovis, hvplot 
#              and matplotlib
#
#  updated 2024/6  matplotlib polycollection using dealiased=False,
#  which makes matplotlib much faster and removes the lines between
#  polygons
#
#  Writing a PNG of size ~ 6400x3200
#  plotting the oRRS18to6 polygons (shaded by area): 
#     holoviews (matplotlib extension)                   26.39s
#     hvplot (matplotlib backend):                       21.38
#     matplotlib's polycollection:                       24.21
#     holoviews (bokeh extension/firefox renderer):      38.37
#     matplotlib/interp to latlon/pcolor                239s
#     holoviews/matplotlib/hv.save (no rasterizations) 2852s
#     holoviews/bokeh/html                             1799s
#
#  plotting ne1024pg2 polygons (Needs ~55GB of memory)
#     holoviews (matplotlib extension)                   103s
#     hvplot (matplotlib backend):                        98s
#     matplotlib's polycollection:                       134s
#     holoviews (bokeh extension/firefox renderer):      122s
#
#
# Notes:
#   hvplot is the fastest.  hvplot is a thin frontend to holoviews, so we
#   should be able to reproduce hvplot speed using holoviews directly.
#   hvplot seems to use a different rasterization approach.  When underresolved
#   (i.e. the 1600x800 image with the oRRS18to6 grid), it has some white artifacts
#   in the plot.  But these go away at 6400x3200
#
#   The other approaches all produce good results at 1600x800 and 6400x3200
#
#   mfa.py currently uses the holoviews/rasterization/matplot extension 
#   approach.
#
    
from netCDF4 import Dataset
import time

import shapely
import pyarrow as pa
import sys
import spatialpandas
import cartopy.crs as ccrs

import numpy as np
import matplotlib
from matplotlib import pyplot
from matplotlib.collections import PolyCollection
from scipy.interpolate import griddata

import holoviews as hv
from holoviews.operation.datashader import rasterize as hds_rasterize
import bokeh

from selenium.webdriver import Firefox
from selenium.webdriver.firefox.options import Options as FirefoxOptions
import traceback

import hvplot.pandas



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
    
    dproj=ccrs.PlateCarree()
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
        xv,yv=np.meshgrid(lon_i,lat_south)
        coords_in  = ccrs.SouthPolarStereo().transform_points(dproj,lon_h,lat_h)
        coords_out = ccrs.SouthPolarStereo().transform_points(dproj,xv.flatten(),yv.flatten())
        data_s = griddata(coords_in[:,0:2], data2d_h, coords_out[:,0:2], method='linear')
        
        data2d_h=data2d[lat>-halo]
        lon_h=lon[lat>-halo]
        lat_h=lat[lat>-halo]
        xv,yv=np.meshgrid(lon_i,lat_north)
        coords_in  = ccrs.NorthPolarStereo().transform_points(dproj,lon_h,lat_h)
        coords_out = ccrs.NorthPolarStereo().transform_points(dproj,xv.flatten(),yv.flatten())
        data_n = griddata(coords_in[:,0:2], data2d_h, coords_out[:,0:2], method='linear')
        
        data_i=np.concatenate((data_s,data_n)).reshape(len(lat_i),len(lon_i))

    elif lat_i[-1]<0:
        # SH only
        data2d_h=data2d[lat<halo]
        lon_h=lon[lat<halo]
        lat_h=lat[lat<halo]
        
        xv,yv=np.meshgrid(lon_i,lat_i)
        coords_in  = ccrs.SouthPolarStereo().transform_points(dproj,lon_h,lat_h)
        coords_out = ccrs.SouthPolarStereo().transform_points(dproj,xv.flatten(),yv.flatten())
        data_i = griddata(coords_in[:,0:2], data2d_h, coords_out[:,0:2], method='linear')

    elif lat_i[0]>0:
        # NH only
        data2d_h=data2d[lat>-halo]
        lon_h=lon[lat>-halo]
        lat_h=lat[lat>-halo]
        
        xv,yv=np.meshgrid(lon_i,lat_i)
        coords_in  = ccrs.NorthPolarStereo().transform_points(dproj,lon_h,lat_h)
        coords_out = ccrs.NorthPolarStereo().transform_points(dproj,xv.flatten(),yv.flatten())
        data_i = griddata(coords_in[:,0:2], data2d_h, coords_out[:,0:2], method='linear')

    else:
        print("Error: interp_to_latlon failed")
        sys.exit(1)
        
    return data_i




def shift_anti_meridian_polygons(lon_poly_coords, lat_poly_coords, eps=40):
    """Shift polygons that are split on the anti-meridian for visualization
    
    Parameters
    ----------
    lon_poly_coords : ndarray(n, v)
        longitudinal coordinates of each vertex of each polygon
    lat_poly_coords : ndarray(n, v)
        latitudinal coordinates of each vertex of each polygon
    eps : float (default 10)
        Tolerance for polygons to shift
    
    Returns
    -------
    polygons : ndarray(n, v, 2)
        Combined longitude and latitude coordinates for all polygons (including shifted ones)
    """
    polygons = np.stack((lon_poly_coords, lat_poly_coords), axis=2)
    diff = np.array(np.max(polygons[:,:,0], axis=1) - np.min(polygons[:,:,0], axis=1) > eps)
    lon_coord_mask = polygons[:,:,0] < eps
    lon_coord_mask[~diff,:] = 0

    polygons[lon_coord_mask,0] = polygons[lon_coord_mask,0] + 360
    return polygons


def polygons_to_geodataframe(lon_poly_coords, lat_poly_coords, data, eps=10):
    """Takes Coordinates for polygons and converts them to a geodataframe for Holoviews Visualization

    Data Structure Transition & Explanation: 
        numpy -> shapely -> pyarrow -> spatialpandas multipolygonarray -> spatialpandas geodataframe

        pyarrow is a Python implementation of Apache Arrow which optimizes data transition from JSON
        to dataframes that gives us the majority of the speedup seen from using the geodataframe 
        instead of just numpy arrays. This function currently only takes milliseconds to run as well.
    
    Parameters
    ----------
    lon_poly_coords : ndarray(n, v)
        longitudinal coordinates of each vertex of each polygon
    lat_poly_coords : ndarray(n, v)
        latitudinal coordinates of each vertex of each polygon
    data : ndarray(n)
        The data to assign to each of the n polygons
    eps : float (default 10)
        Tolerance for polygons to shift
    
    Returns
    -------
    polygons : spatialpandas.GeoDataFrame columns: {'geometry': polygons, 'faces': data}
        Geodataframe of all cell polygons with data linked to each cell
    """

    polygons = shift_anti_meridian_polygons(lon_poly_coords, lat_poly_coords)
    geo = shapely.polygons(polygons)

    arr_flat, part_indices = shapely.get_parts(geo, return_index=True)
    offsets1 = np.insert(np.bincount(part_indices).cumsum(), 0, 0)
    arr_flat2, ring_indices = shapely.get_rings(arr_flat, return_index=True)
    offsets2 = np.insert(np.bincount(ring_indices).cumsum(), 0, 0)
    coords, indices = shapely.get_coordinates(arr_flat2, return_index=True)
    offsets3 = np.insert(np.bincount(indices).cumsum(), 0, 0)
    coords_flat = coords.ravel()
    offsets3 *= 2

    _parr3 = pa.ListArray.from_arrays(pa.array(offsets3), pa.array(coords_flat))
    _parr2 = pa.ListArray.from_arrays(pa.array(offsets2), _parr3)
    parr = pa.ListArray.from_arrays(pa.array(offsets1), _parr2)

    polygons = spatialpandas.geometry.MultiPolygonArray(parr)
    gdf = spatialpandas.GeoDataFrame({'geometry': polygons})
    gdf = gdf.assign(faces = data)
    return gdf







# read in polygons from scrip file:
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
#name="/Users/mataylo/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne256pg2.scrip.nc"
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne1024pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne256pg2.scrip.nc"
name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne1024pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/ocean.oRRS18to6v3.scrip.181106.nc"
#name="/Users/mt/scratch1/mapping/grids/ocean.oRRS18to6v3.scrip.181106.nc"
#name="/Users/mataylo/scratch1/mapping/grids/ocean.oRRS18to6v3.scrip.181106.nc"


file1 = Dataset(name,"r")
ncols = file1.dimensions["grid_size"].size
print("ncols = ",ncols)
clat1 = file1.variables["grid_center_lat"][:]
clon1 = file1.variables["grid_center_lon"][:]
xlat = file1.variables["grid_corner_lat"][:,:]
xlon = file1.variables["grid_corner_lon"][:,:]
area = file1.variables["grid_area"][:]
if "radian" in file1.variables["grid_corner_lat"].units.lower():
    xlat *= 180/np.pi
    xlon *= 180/np.pi



mn=float(min(area))
mx=float(max(area))
clev=(mn,mx)
colormap='Spectral'
print(f"poly_plot(): plotting {len(area)} cells. data min/max= {mn:.3},{mx:.3}")

# center plot at lon=0,lat=0:
proj=ccrs.PlateCarree()
xpoly  = proj.transform_points(proj, xlon, xlat)

# convert polygons go geodataframe (for holoviews routines)
xlon=xpoly[:,:,0]
xlat=xpoly[:,:,1]


# holovies resolution:
height=800
width=1600
# mpl resolution:
dpi=300    # only used by MPL direct code. dpi=300 results in  1548x804 image

# bump up the resolution:
height=4*height
width=4*width
dpi=4*dpi

# global plot
xlim=(-180.,180.)
ylim=(-90.,90.)





################################################################################
#
# holoviews/matplotlib/hds_rasterize
#
#
################################################################################
print("holoviews/matplotlib/hds_rasterize... ",end='')
start= time.time()
hv.extension('matplotlib')
gdf = polygons_to_geodataframe(np.ma.getdata(xlon),np.ma.getdata(xlat), np.ma.getdata(area))
hv_polys = hv.Polygons(gdf, vdims=['faces'])
hv_polys.opts(color='faces')
hv_polys.opts(colorbar=False)

# cant figure out how to change dpi, seems fixed at 72
r= hds_rasterize(hv_polys,width=width,height=height)  
r.opts(xlim=(-180.,180))
r.opts(ylim=(-90.,90))
r.opts(clim=clev)
r.opts(cmap=colormap)
r.opts(data_aspect=1)
r.opts(xaxis=None, yaxis=None)
r.opts(fig_inches=width/72)   # this works for matplotlib backend



fig=hv.render(r)     # adding dpi here doesn't change anything, width not accepted
fig.savefig("hv-mpl-rasterize.png", bbox_inches='tight')
end=time.time()
print(f"{end-start:.2f}s")




################################################################################
#
# hvplot (wrapper to holowviews)
#
################################################################################
print("hvplot/matplotlib/savefig... ",end="")
start= time.time()
hv.extension('matplotlib') 

gdf = polygons_to_geodataframe(np.ma.getdata(xlon),np.ma.getdata(xlat), np.ma.getdata(area))
plot2 = gdf.hvplot.polygons(cmap=colormap,clim=clev,
              color='faces',
           width=width, height=height, xlim=xlim, ylim=ylim, 
           edgecolor='none',alpha=1,
                            xaxis=None, yaxis=None,
             data_aspect=1, colorbar=False, rasterize=True)
# hv seems to work with none, dont need  edgecolor='face'?

fig = hv.render(plot2)
fig.savefig("hvplot-mpl.png",  bbox_inches='tight')
end=time.time()
print(f"{end-start:.2f}s")




################################################################################
#
#matplotlib's polycollection:
#
################################################################################
print("matplotlib/polycollection... ",end='')
start= time.time()

# adjust cells into polycollection format:
xpoly = shift_anti_meridian_polygons(xpoly[:,:,0],xpoly[:,:,1])

ax = pyplot.axes(projection=ccrs.PlateCarree())
ax.set_global()

# antialized=False to remove white lines between cells.
p = matplotlib.collections.PolyCollection(xpoly, array=area,
    edgecolor='none',alpha=1,antialiased=False)

p.set_clim(clev)
p.set_cmap(colormap)
ax.add_collection(p)

matplotlib.pyplot.savefig('mpl-pc.png',dpi=dpi,orientation="portrait",bbox_inches='tight')
end= time.time()
print(f"{end-start:.2f}s")

################################################################################
#
# holoviews/bokeh/firefox 
#
################################################################################
print("holoviews/bokeh/Firefox... ",end='')
start= time.time()
hv.extension('bokeh') # need to load extension before setting options
gdf = polygons_to_geodataframe(np.ma.getdata(xlon),np.ma.getdata(xlat), np.ma.getdata(area))
hv_polys = hv.Polygons(gdf, vdims=['faces']).opts(color='faces')
rasterized = hds_rasterize(hv_polys,height=height, width=width)

out_plot = rasterized.opts(data_aspect=1,  height=height, width=width,
                               xlim=xlim, ylim=ylim, \
                               cmap=colormap,clim=clev,\
                              xaxis=None, yaxis=None,
                           colorbar=False)


        
web_driver = None
options = FirefoxOptions()
options.add_argument('--headless') # Argument for web driver with no user interface

# Web drivers become orphan processes if quit() is not called
try:
    web_driver = Firefox(options=options)
    fig = hv.render(out_plot)
    bokeh.io.export_png(fig, filename="holoviews-bokey-firefox.png", height=height, width=width, webdriver=web_driver)
    web_driver.quit()
    
except Exception:
    traceback.print_exception(*sys.exc_info())
    if web_driver is not None:
        web_driver.quit()
        sys.exit()

end=time.time()
print(f"{end-start:.2f}s")



sys.exit(0)
#
#the approaches below take 10x longer than the above
#probably related to when the polygons are rasterized
#

################################################################################
#
# inteprolate to lat/lon
#
#
################################################################################
print("matplotlib/interpolate to lat/lon... ",end='')
start= time.time()
nlat=height // 2
nlon=2*nlat
lon_i = np.linspace(0, 360, nlon,endpoint=False)  
dlat2=90./nlat
lat_i = np.linspace(-90+dlat2, 90-dlat2, nlat)  
data_i=interp_to_latlon(area,clat1,clon1,lat_i,lon_i)

ax = pyplot.axes(projection=ccrs.PlateCarree())
ax.set_global()
pl=ax.pcolormesh(lon_i, lat_i, data_i,vmin=mn,vmax=mx, cmap=colormap)
pyplot.savefig('mpl-latlon.png',dpi=dpi,orientation="portrait",bbox_inches='tight')

end=time.time()
print(f"{end-start:.2f}s")





################################################################################
#
# holoviews/matplotlib
#
# terribly slow:
#
################################################################################
print("holoviews/matplotlib/hv.save... ",end='')
start= time.time()
hv.extension('matplotlib') 
gdf = polygons_to_geodataframe(np.ma.getdata(xlon),np.ma.getdata(xlat), np.ma.getdata(area))
hv_polys = hv.Polygons(gdf, vdims=['faces'])
hv_polys.opts(color='faces')
hv_polys.opts(cmap=colormap)
hv_polys.opts(data_aspect=1)
hv_polys.opts(edgecolor='none')    
hv_polys.opts(alpha=1)
hv_polys.opts(xlim=(-180.,180))
hv_polys.opts(ylim=(-90.,90))
hv_polys.opts(clim=clev)
#hv_polys.opts(width=width,height=height)  # "not valid for matplotlib backend
hv_polys.opts(colorbar=False)
#hv_polys.opts(xlabel='', ylabel='')
#hv_polys.opts(fontscale=0.25)
#hv_polys.opts(fontsize={'xticks': 0,'yticks': 0})
hv_polys.opts(xaxis=None, yaxis=None)
#hv_polys.opts(rasterize=True) # not valid


# tweak dpi to get close to 1600x800 plot
hv.save(hv_polys, filename="hv-mpl-save.png") 
end=time.time()
print(f"{end-start:.2f}s")

# this is always about 2x slower than hv.save
#print("holoviews/matplotlib/savefig:")
#start= time.time()
#fig1 = hv.render(hv_polys)  # render via matplotlib
#fig1.savefig("area-mpl-savefig.png", dpi=dpi, bbox_inches='tight')
#end=time.time()
#print(f"holoviews/matplotlib/savefig: {end-start:.2f} s")



################################################################################
#
# holoviews to html   (the only scalable vector graphics version)
#
################################################################################
print("calling holoviews/bokeh/html... ",end='')
start= time.time()
xlim=(-180.,180.)
ylim=(-90.,90.)
hv.extension('bokeh') # need to load extension before setting options
gdf = polygons_to_geodataframe(np.ma.getdata(xlon),np.ma.getdata(xlat), np.ma.getdata(area))
hv_polys = hv.Polygons(gdf, vdims=['faces'])

hv_polys.opts(colorbar=False,cmap=colormap,clim=clev,
              color='faces',
              xlim=xlim,ylim=ylim,data_aspect=1,
                      line_width=0,line_alpha=1,alpha=1,
              width=width,height=height)
hv.save(hv_polys, 'hv-bokeh-html.html')

end=time.time()
print(f"{end-start:.2f}s")

# requires firefox and geckodriver or chormium and chormedriver available in PATH:
#start= time.time()
#hv.save(hv_polys, 'area-hv-bokeh-png.png')
#end=time.time()
#print(f"holoviews/bokeh/hv.save(png): {end-start:.2f} s")


