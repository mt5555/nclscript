import numpy as np
#
import time
import sys
from netCDF4 import Dataset

from plotpoly_mpl import plotpoly
from plotpoly_hv import plotpoly as plotpoly_hv
import cartopy.crs as ccrs

# read in polygons from scrip file:
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
name="/Users/mataylo/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
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
Rearth_km = 6378.1                # radius of earth, in km
res=Rearth_km*np.sqrt(area)

alpha=np.ones_like(res)
alpha[res>150]=.5
alpha[res>160]=0

#proj=ccrs.PlateCarree()
proj=ccrs.Robinson()
#clat=40; clon=-60;  proj = ccrs.Orthographic(central_latitude=clat, central_longitude=clon) 
#print(proj.srs)


#plotpoly(xlat,xlon,res,"res.png",title="resolution (km)",alpha=alpha)
plotpoly_hv(xlat,xlon,res,"res.png",title="resolution (km)",proj=proj,alpha=alpha)


