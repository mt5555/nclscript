import numpy as np
#
# hv-example.py
#       python code using holoviews to make an unstructured grid plot
#       plots cell center data, coloring each cell based on data values
#       cells are represented as polygons, read in from grid "scrip" file
#
#       for simplicity, this script plots cell area, since this data is
#       contained in the scrip file.
#
#
import time
import sys
from netCDF4 import Dataset

import shapely
import pyarrow as pa
import spatialpandas
import cartopy.crs as ccrs
import matplotlib

import holoviews as hv
from holoviews.operation.datashader import rasterize as hds_rasterize



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
name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne256pg2.scrip.nc"
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne1024pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne256pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne1024pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/ocean.oRRS18to6v3.scrip.181106.nc"
#name="/Users/mt/scratch1/mapping/grids/ocean.oRRS18to6v3.scrip.181106.nc"


file1 = Dataset(name,"r")
ncols = file1.dimensions["grid_size"].size
print("ncols = ",ncols)
#clat1 = file1.variables["grid_center_lat"][:]
#clon1 = file1.variables["grid_center_lon"][:]
xlat = file1.variables["grid_corner_lat"][:,:]
xlon = file1.variables["grid_corner_lon"][:,:]
area = file1.variables["grid_area"][:]
if "adian" in file1.variables["grid_corner_lat"].units:
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
fig.savefig("hv-example.png", bbox_inches='tight')




