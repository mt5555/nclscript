import numpy as np
#
import time
import sys
from netCDF4 import Dataset

import matplotlib
import matplotlib.pylab
from plotpoly_mpl import plotpoly
from plotpoly_hv import plotpoly as plotpoly_hv
import cartopy.crs as ccrs

# read in polygons from scrip file:
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
#name="/Users/mataylo/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne120pg2.scrip.nc"
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne256pg2.scrip.nc"
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne1024pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne256pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne1024pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/ocean.oRRS18to6v3.scrip.181106.nc"
#name="/Users/mt/scratch1/mapping/grids/ocean.oRRS18to6v3.scrip.181106.nc"


file1 = Dataset(name,"r")
ncols = file1.dimensions["grid_size"].size
#clat1 = file1.variables["grid_center_lat"][:]
#clon1 = file1.variables["grid_center_lon"][:]
xlat = file1.variables["grid_corner_lat"][:,:]
xlon = file1.variables["grid_corner_lon"][:,:]
area = file1.variables["grid_area"][:]
Rearth_km = 6378.1                # radius of earth, in km
res=Rearth_km*np.sqrt(area)


#proj=ccrs.PlateCarree()
#proj=ccrs.Robinson()
clat=40; clon=-60;  proj = ccrs.Orthographic(central_latitude=clat, central_longitude=clon) 
print(proj.srs)

# adding alpha array with data:
#  matplotlib:  works, but doesn't apply to cell edges
#  hv:  ignored by rasterize
#alpha=np.ones_like(res)
#alpha[res>150]=.5
#alpha[res>160]=0
alpha=1

# adding alpha to colormap:
#  matplotlib:  works, but doesn't apply to cell edges
#  hv:  works!
cmap=matplotlib.pylab.cm.plasma
my_cmap = cmap(np.arange(cmap.N))
#my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
my_cmap[0:round(cmap.N/4),-1] = 0
my_cmap[round(cmap.N/4):round(cmap.N/2),-1]=0.5
my_cmap = matplotlib.colors.ListedColormap(my_cmap)
colormap=my_cmap

t1= time.time()
plotpoly(xlat,xlon,res,"resmpl.png",title="resolution (km)",proj=proj,colormap=colormap)
t2= time.time()
print(f"mpl polycollection: {t2-t1:.2f}s")
plotpoly_hv(xlat,xlon,res,"reshv.png",title="resolution (km)",proj=proj,colormap=colormap)
t3= time.time()
print(f"hv polycollection:  {t3-t2:.2f}s")

