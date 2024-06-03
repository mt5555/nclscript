import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

# Load the data from the file
data = nc.Dataset('/ascldap/users/ljmaloy/data/TEMPEST_ne30pg2.scrip.nc')

# Extract the longitude and latitude arrays
grid_center_lon = data.variables['grid_center_lon'][:]
grid_center_lat = data.variables['grid_center_lat'][:]
grid_corner_lon = data.variables['grid_corner_lon'][:]
grid_corner_lat = data.variables['grid_corner_lat'][:]

def shift_anti_meridian_polygons(lon_poly_coords, lat_poly_coords, eps=40):
    polygons = np.stack((lon_poly_coords, lat_poly_coords), axis=2)
    diff = np.array(np.max(polygons[:,:,0], axis=1) - np.min(polygons[:,:,0], axis=1) > eps)
    lon_coord_mask = polygons[:,:,0] < eps
    lon_coord_mask[~diff,:] = 0
    polygons[lon_coord_mask,0] += 360
    return polygons

# Choose projection
#proj = ccrs.PlateCarree()
#proj = ccrs.Orthographic(central_latitude=40, central_longitude=-30)
#proj = ccrs.Mollweide()
proj = ccrs.Robinson()
#proj = ccrs.Stereographic(central_latitude=90, central_longitude=-45)

# Plot 1: Using lat/lon data directly
fig1, ax1 = plt.subplots(subplot_kw={'projection': proj})
ax1.coastlines(resolution='110m')

# Plotting cell centers
ax1.scatter(grid_center_lon, grid_center_lat, s=1, transform=ccrs.PlateCarree(), color='blue')

# Plotting polygons
if proj == ccrs.PlateCarree():
    transformed_polygons = shift_anti_meridian_polygons(grid_corner_lon, grid_corner_lat)
else:
    transformed_polygons = np.stack((grid_corner_lon, grid_corner_lat), axis=2)

for i in range(len(grid_corner_lon)):
    ax1.add_patch(plt.Polygon(transformed_polygons[i], fill=None, edgecolor='black', linewidth=0.5, transform=ccrs.PlateCarree()))
ax1.set_title('Direct Lat/Lon Plot')

# Plot 2: Transforming data before plotting
fig2, ax2 = plt.subplots(subplot_kw={'projection': proj})
ax2.coastlines(resolution='110m')

# Transforming cell centers
transformed_points = proj.transform_points(ccrs.PlateCarree(), grid_center_lon, grid_center_lat)
ax2.scatter(transformed_points[:, 0], transformed_points[:, 1], s=1, color='red')

# Transforming and plotting polygons
for i in range(len(grid_corner_lon)):
    transformed_poly = proj.transform_points(ccrs.PlateCarree(), grid_corner_lon[i], grid_corner_lat[i])
    ax2.add_patch(plt.Polygon(transformed_poly[:, :2], fill=None, edgecolor='purple', linewidth=0.5))
ax2.set_title('Transformed Data Plot')

plt.show()
input("Press a key to continue...")
