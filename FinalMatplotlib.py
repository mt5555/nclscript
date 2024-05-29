import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
 
# Load the data from the file
data = nc.Dataset('TEMPEST_ne30pg2.scrip.nc')
 
# Extract the longitude and latitude arrays
grid_center_lon = data.variables['grid_center_lon'][:]
grid_center_lat = data.variables['grid_center_lat'][:]
grid_corner_lon = data.variables['grid_corner_lon'][:]
grid_corner_lat = data.variables['grid_corner_lat'][:]
 
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
 
# Uncomment for different projections
#proj = ccrs.Orthographic(central_latitude=40, central_longitude=-30)
proj = proj=ccrs.PlateCarree()
#proj = ccrs.Mollweide()
#proj = ccrs.Robinson()
#proj = ccrs.Stereographic(central_latitude=90, central_longitude=-45)
 
if isinstance(proj, ccrs.PlateCarree):
    xpoly2 = proj.transform_points(proj, grid_corner_lon, grid_corner_lat)
    xlon = xpoly2[:,:,0]
    xlat = xpoly2[:,:,1]
    xpoly2 = shift_anti_meridian_polygons(xlon, xlat)
 
# Create the plot with Orthographic projection
fig, ax = plt.subplots(subplot_kw={'projection': proj})
 
# Add continental outlines for Orthographic
ax.coastlines(resolution='110m')
 
# Plot the cells
for i in range(len(grid_center_lon)):
    if isinstance(proj, ccrs.PlateCarree):
        cell = list(zip(xpoly2[i][:, 0], xpoly2[i][:, 1]))
    else:
    cell = list(zip(grid_corner_lon[i], grid_corner_lat[i]))
    ax.add_patch(plt.Polygon(cell, fill=None, edgecolor='black', linewidth=0.5, transform=ccrs.PlateCarree()))
 
# Transform the data
transformed_lon, transformed_lat, transformed_polygon = proj.transform_points(ccrs.PlateCarree(), grid_center_lon, grid_center_lat).T
 
# Scatter plot of transformed grid points
ax.scatter(transformed_lon, transformed_lat, s=1, transform=proj)
 
# Set labels and title
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Transformed Grid Points')
 
# Display the plot
plt.show()
 
input("Press a key to continue...")
