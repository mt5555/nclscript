import numpy as np

import shapely
import pyarrow as pa
import sys
import spatialpandas
import cartopy.crs as ccrs
import cartopy.feature as cf

from PIL import Image   # needed to load JPG background image

import geoviews as gv
import geoviews.feature as gf
from holoviews.operation.datashader import datashade,rasterize 
import time

def shift_anti_meridian_polygons(polygons, eps=40):
    #shift polygons that are split on the anti-meridian for visualization

    diff = np.array(np.max(polygons[:,:,0], axis=1) - np.min(polygons[:,:,0], axis=1) > eps)
    lon_coord_mask = polygons[:,:,0] < eps   # all polygons on left edge
    lon_coord_mask[~diff,:] = 0              # mask=0 for subset of left polygons which are not cut
    polygons_new=polygons[diff,:,:]            # set of all split polygons
    polygons[lon_coord_mask,0] = polygons[lon_coord_mask,0] + 360

    lon_coord_mask = polygons_new[:,:,0] > eps  # coords on right side
    polygons_new[lon_coord_mask,0] = polygons_new[lon_coord_mask,0] - 360
    # also return polygons_new, and "diff", so we can extract
    # data_new = data[diff] 
    return [polygons,polygons_new,diff]


def polygons_to_geodataframe(polygons, data, eps=10):
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



def plotpoly(lat_poly_coords, lon_poly_coords, data, filepath=None, title='',
              proj=ccrs.PlateCarree(),width=4000, height=1800, 
              xlim=(-180.,180), ylim=(-90.,90.),
              clim=None,colormap=None,mask=1,alpha=1
):
    if lon_poly_coords.shape != lat_poly_coords.shape:
        print(f"Dimension mismatch between longitude: {lon_poly_coords.shape}, and latitude: {lat_poly_coords.shape}")
        return
    elif len(lon_poly_coords) != len(data):
        print(f"Dimension mismatch between number of cells: {len(lon_poly_coords)} and number of data points: {len(data)}")
        return

    # if mask present, remove masked cells
    if not np.isscalar(mask):
        data=data[mask]
        lon_poly_coords = lon_poly_coords[mask,:]
        lat_poly_coords = lat_poly_coords[mask,:]
        #count = sum(1 for x in mask if x)

    # convert to degrees, if necessary
    if np.max(np.abs(lat_poly_coords))<1.1*np.pi:
        lat_poly_coords=lat_poly_coords*180/np.pi
        lon_poly_coords=lon_poly_coords*180/np.pi
        
    mn=float(min(data))
    mx=float(max(data))
    print(f"plotpoly(): {len(data)} cells. data min/max= {mn:.3},{mx:.3} {title}")
    if clim == None:
        clim=(mn,mx)
    print(f"clim = ({clim[0]:.3},{clim[1]:.3})")
    if colormap==None:
        if mn*mx < 0: colormap='Spectral'
        else: colormap='plasma'


    # transform into desired coordinate system:
    xpoly  = proj.transform_points(ccrs.PlateCarree(), lon_poly_coords, lat_poly_coords)
    #xpoly  = gv.operation.project(projection=proj, corners)

    # fix and duplicate cut cells
    if "proj=eqc" in proj.srs:
        # duplicate cut polygons on left and right edge of plot
        [xpoly,xpoly_new,mask_new] = shift_anti_meridian_polygons(xpoly[:,:,0:2])
        corners=np.concatenate((xpoly[:,:,0:2],xpoly_new[:,:,0:2]),axis=0)
        data=np.concatenate((data,data[mask_new]),axis=0)
        if not np.isscalar(alpha):
            alpha=np.concatenate((alpha,alpha[mask_new]),axis=0)
    if "proj=robin" in proj.srs:
        # remove all cut polygons
        eps=40*1e5
        mask_keep = np.array(np.max(xpoly[:,:,0], axis=1) - np.min(xpoly[:,:,0], axis=1) < eps)
        corners=xpoly[mask_keep,:,0:2]
        data=data[mask_keep]
        if not np.isscalar(alpha):
            alpha=alpha[mask_keep]
    if "proj=ortho" in proj.srs:
        #remove non-visible points:
        mask_keep =  np.all(np.isfinite(xpoly),axis=(1,2))
        corners = xpoly[mask_keep,:,0:2]
        data=data[mask_keep]
        if not np.isscalar(alpha):
            alpha=alpha[mask_keep]

    gdf = polygons_to_geodataframe(np.ma.getdata(corners[:,:,:]), np.ma.getdata(data[:]))

    cbar_opts={}
    #cbar_opts={'width': round(.02*width)}
    #cbar_opts={'label': "km"}

    gv.extension('matplotlib') # need to load extension before setting options
    hv_polys = gv.Polygons(gdf, vdims=['faces'],crs=proj).opts(color='faces')
    hv_polys.opts(projection=proj, global_extent=True)
    #r.opts(xlim=(-180.,180))
    #r.opts(ylim=(-90.,90))


    rasterized = rasterize(hv_polys,height=height, width=width)
    rasterized.opts(cmap=colormap,colorbar=True,colorbar_opts=cbar_opts)
    rasterized.opts(clim=clim)
    rasterized.opts(fontscale=10)
    rasterized.opts(title=title)
    rasterized.opts(xlabel='', ylabel='', clabel='')


    #background = gv.Overlay([gf.coastline,gf.ocean,gf.land])
    #r = background * rasterized
    

    #https://www.color-hex.com/color-palette/1021516  ocean: #004589  land: #72601b 
    # blue from 2012 ANL video:  rgb(48 62 141) = #303E8D   nice blue
    # blue from me:  #01013f                                dark blue
    # blue from world iamge      rgb(2  5  20) =#020514     almost black
    # background = gf.ocean.opts(facecolor='#01013f') * gf.land.opts(facecolor='#958258') 

    #background: image
    Image.MAX_IMAGE_PIXELS = 233280000
    image_path = 'world.topo.200408.3x5400x2700.png'
    #image_path = 'world.topo.200408.3x21600x10800.png'
    print(f"background={image_path}")
    img_data = np.flipud( np.array(Image.open(image_path)) )
    bounds = (-180, -90, 180, 90)  # Assuming the image covers the whole globe
    background=gv.RGB((np.linspace(-180, 180, img_data.shape[1]),
                 np.linspace(-90, 90, img_data.shape[0]),
                 img_data[..., 0], img_data[..., 1], img_data[..., 2]),
                  bounds=bounds,crs=ccrs.PlateCarree()).opts(projection=proj)
 
    
    r = background * rasterized    
    r.opts(fig_inches=width/72)   
    r.opts(data_aspect=1)
    #r.opts(fig_size=100)  # scaling factor

    print("gv rendering...")
    fig=gv.render(r)
    
    if (filepath!=None):
        print(f"writing: {filepath}")
        fig.savefig(filepath, bbox_inches='tight')

