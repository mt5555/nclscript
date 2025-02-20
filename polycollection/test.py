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
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne120pg2.scrip.nc"
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne256pg2.scrip.nc"
#name="/Users/mt/scratch1/mapping/grids/TEMPEST_ne1024pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne256pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/TEMPEST_ne1024pg2.scrip.nc"
#name="/ascldap/users/mataylo/scratch1/mapping/grids/ocean.oRRS18to6v3.scrip.181106.nc"
#name="/Users/mt/scratch1/mapping/grids/ocean.oRRS18to6v3.scrip.181106.nc"

#name="/home/ac.mtaylor/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
#dname="/home/ac.mtaylor/scratch1/viz/output.scream.decadal.6hourlyINST_ne30pg2.INSTANT.nhours_x6.1994-10-26-21600.nc"

name="/home/ac.mtaylor/scratch1/mapping/grids/TEMPEST_ne1024pg2.scrip.nc"
dname="/home/ac.mtaylor/scratch1/viz/output.scream.decadal.1hourlyINST_ne1024pg2.INSTANT.nhours_x1.1994-11-15-03600.nc"



print(f"input: {dname}")
file1 = Dataset(name,"r")
ncols = file1.dimensions["grid_size"].size
#clat1 = file1.variables["grid_center_lat"][:]
#clon1 = file1.variables["grid_center_lon"][:]
xlat = file1.variables["grid_corner_lat"][:,:]
xlon = file1.variables["grid_corner_lon"][:,:]

#area = file1.variables["grid_area"][:]
#Rearth_km = 6378.1                # radius of earth, in km
#res=Rearth_km*np.sqrt(area)

datafile=Dataset(dname,"r")
#res = datafile.variables["qv_2m"][0,:]  # first snapshot in file
res = datafile.variables["VapWaterPath"][0,:]  # first snapshot in file


#proj=ccrs.PlateCarree()
proj=ccrs.Robinson()
#clat=40; clon=-60;  proj = ccrs.Orthographic(central_latitude=clat, central_longitude=clon) 
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
# cmap=matplotlib.pylab.cm.plasma
# my_cmap = cmap(np.arange(cmap.N))
# # add transparency:
# my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
# my_cmap[0:round(cmap.N/4),-1] = 0
# my_cmap[round(cmap.N/4):round(cmap.N/2),-1]=0.5
# my_cmap = matplotlib.colors.ListedColormap(my_cmap)

# Create the custom colormap
#cmap = matplotlib.pyplot.get_cmap("Greys_r")(np.linspace(0, 1, 256))
#cmap = cmap[40:]
#cmap[:int(len(cmap)/4), 3] = np.linspace(0, 1, int(len(cmap)/4))
#my_cmap = matplotlib.colors.ListedColormap(cmap)


# 100% white, with transparancy
#cmap = matplotlib.pyplot.get_cmap("Greys_r")(np.linspace(0, 1, 256))
#cmap[:,0]=1 ; cmap[:,1]=1 ; cmap[:,2]=1
#cmap[:, 3] = np.linspace(0, 1, len(cmap))**(1.0)
#my_cmap = matplotlib.colors.ListedColormap(cmap)


cmap = matplotlib.pyplot.get_cmap("RdYlBu")(np.linspace(0, 1, 256))
cmap[0:32, 3] = 0
cmap[32:96,3] = np.linspace(0, 1, 64)
my_cmap = matplotlib.colors.ListedColormap(cmap)



#t1= time.time()
#plotpoly(xlat,xlon,res,"resmpl.png",title="resolution (km)",proj=proj,colormap=my_cmap)
t2= time.time()
#print(f"mpl polycollection: {t2-t1:.2f}s")
plotpoly_hv(xlat,xlon,res,"reshv.png",clim=(-10.,70.),title="resolution (km)",proj=proj,colormap=my_cmap)
t3= time.time()
print(f"hv polycollection:  {t3-t2:.2f}s")



#
#  hv with world 5400x2700 background and ne30 data:  33s
#  hv with world 21600x10800 background and ne30 data: 135s
#
##  hv with world 5400x2700 background and ne1024 data:  200s   max memory:  ~30GB
#
