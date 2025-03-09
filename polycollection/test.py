import numpy as np
#
import time
import sys
from netCDF4 import Dataset
from pathlib import Path
import pickle


import matplotlib
import matplotlib.pylab
import geoviews as gv
import geoviews.feature as gf
from plotpoly_mpl import plotpoly
from plotpoly_hv import plotpoly as plotpoly_hv
import cartopy.crs as ccrs
from PIL import Image   # needed to load JPG background image

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
#dname="/home/ac.mtaylor/scratch1/viz/1994/output.scream.decadal.6hourlyINST_ne30pg2.INSTANT.nhours_x6.1994-10-26-21600.nc"

name="/home/ac.mtaylor/scratch1/mapping/grids/TEMPEST_ne1024pg2.scrip.nc"
dname="/home/ac.mtaylor/scratch1/viz/1995-testb2/output.scream.decadal.1hourlyINST_ne1024pg2.INSTANT.nhours_x1.1995-08-20-03600.nc"



print(f"input: {dname}")
file1 = Dataset(name,"r")
ncols = file1.dimensions["grid_size"].size
#clat1 = file1.variables["grid_center_lat"][:]
#clon1 = file1.variables["grid_center_lon"][:]
xlat = file1.variables["grid_corner_lat"][:,:]
xlon = file1.variables["grid_corner_lon"][:,:]

#area = file1.variables["grid_area"][:]
#Rearth_km = 6378.1                # radius of earth, in km
#var=Rearth_km*np.sqrt(area)

datafile=Dataset(dname,"r")
#vas = datafile.variables["qv_2m"][0,:]  # first snapshot in file

varname="Preciptable Water"   ;  varnamef="VapWaterPath"  ; pngname='tmq'
#varname="LW"   ;  varnamef="LW_flux_up_at_model_top"      ; pngname='lw'
#varname="SW"   ;  varnamef="SW_flux_up_at_model_top"       ; pngname='sw'
#varname="Vapor (2m)"   ;  varnamef="qv_2m"                ; pngname='qv2m'
#varname="Precip"   ;  varnamef="precip_total_surf_mass_flux"  ; pngname='prec'

idx=18
var = datafile.variables[varnamef][idx,:]  # last snapshot in file
dtime = datafile.variables["time"][idx]  # last snapshot in file
print("time=",dtime)
pngname = f"{pngname}-{dtime:.2f}"
#print("nan_count =",np.count_nonzero(np.isnan(var)))

pn=3
if pn==1:
    proj=ccrs.PlateCarree() ; projname="latlon0"
    wres=2000 ; hres=round(wres/2)   # ne1024/latlon speckling at wres=10k  uggh
    dpi=1200                         # used by matplotlib version
if pn==2:
    proj=ccrs.Robinson()   ; projname="robinson0"
    wres=10000 ; hres=round(wres/2)
    dpi=1200
if pn==3:
    clat=30.; clon=-60.;
    clon = 180. - np.mod(dtime,1)*360.   # follow the sun
    proj = ccrs.Orthographic(central_latitude=clat, central_longitude=clon)
    projname=f"ortho_{clat:.0f}_{clon:.0f}"
    wres=2000 ; hres=wres    # ne1024/ortho needs wres>2000 to avoid speckling
    dpi=1200

print(proj.srs)




# adding alpha array with data:
#  matplotlib:  works, but doesn't apply to cell edges
#  hv:  ignored by rasterize
#alpha=np.ones_like(var)
#alpha[var>150]=.5
#alpha[var>160]=0

 



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

# create 3 colormaps:
#    my_cmap         =  standard, without transparency
#    my_cmap_alpha   =  colormap with alpha channel
#    my_cmap_mask    =  colormap that will produce an image with alpha values, for a mask
#
if varname=="SW":
    clim=(-1.,900.)
    cmap = matplotlib.pyplot.get_cmap("Greys_r")(np.linspace(0, 1, 256))
    my_cmap = matplotlib.colors.ListedColormap(cmap)

    # 100% white, with transparancy
    cmap=cmap.copy()  #  make a new copy otherwise we change my_cmap
    alpha=np.linspace(0, 1, len(cmap))
    cmap[:,0]=1 ; cmap[:,1]=1 ; cmap[:,2]=1
    cmap[:, 3] = alpha
    #cmap[0,0:3]=0  # black
    #cmap[0,3:4]=0.70  # dark
    my_cmap_alpha = matplotlib.colors.ListedColormap(cmap)


# GOOD SETTINGS FOR TMQ:
if varname=="Preciptable Water":
    clim=(0.,90.)
    cmap = matplotlib.pyplot.get_cmap("RdBu_r")(np.linspace(0.1, .95, 256))
    my_cmap = matplotlib.colors.ListedColormap(cmap)   # make sure to use values, not reference

    cmap=cmap.copy()  #  make a new copy otherwise we change my_cmap
    opaque=round(.55*my_cmap.N)
    cmap[0:opaque, 3] = np.linspace(0, 1, opaque)
    my_cmap_alpha = matplotlib.colors.ListedColormap(cmap)
    alpha=cmap[:,3]



# create a colormap with colors = alpha 
cmap=cmap.copy()  #  make a new copy otherwise we change my_cmap
cmap[:,0]=alpha
cmap[:,1]=alpha
cmap[:,2]=alpha
cmap[:,3]=1
my_cmap_mask = matplotlib.colors.ListedColormap(cmap,'mt3')



#
#  MPL plot
#
if True:
    t1= time.time()
    plotpoly(xlat,xlon,var,pngname,title=varname,proj=proj,colormap=my_cmap,colormap_mask=my_cmap_mask,clim=clim,dpi=dpi)
    t2= time.time()
    print(f"mpl polycollection: {t2-t1:.2f}s")
    exit(0)


#
#  HV plot
#
#
# for holoviews, convert background image into new projection
# takes about 20min, so save results for next time
#
gv.extension('matplotlib') # need to load extension before setting options
bg=1
if bg==1:
    projname=f"world-hr-{projname}.pickle"
    p=Path(projname)
    if p.exists() and p.is_file():
        print(f"loading background image {p}")
        with open(projname, 'rb') as file:
            background = pickle.load(file)
    else:
        Image.MAX_IMAGE_PIXELS = 233280000
        # image_path = 'world.topo.200408.3x5400x2700.png'
        image_path = 'world.topo.200408.3x21600x10800.png'
        img_data = np.array(Image.open(image_path)) 
        bounds = (-180, 90, 180, -90)  # Assuming the image covers the whole globe
        background=gv.RGB((np.linspace(-180, 180, img_data.shape[1]),
                           np.linspace(90, -90, img_data.shape[0]),
                           img_data[..., 0], img_data[..., 1], img_data[..., 2]),
                          bounds=bounds,crs=ccrs.PlateCarree())  #.opts(projection=proj)
        if pn==1:
            # for latlon projections, dont need to project:
            background.opts(projection=proj)
        else:
            print(f"projecting background image {p}")
            background = gv.project(background,projection=proj)
            print("saving projected background")
            with open(projname, 'wb') as file:
                pickle.dump(background, file)
if bg==2:
    background = gv.Overlay([gf.coastline,gf.ocean,gf.land])
if bg==3:
    #https://www.color-hex.com/color-palette/1021516  ocean: #004589  land: #72601b 
    # blue from 2012 ANL video:  rgb(48 62 141) = #303E8D   nice blue
    # blue from me:  #01013f                                dark blue
    # blue from world iamge      rgb(2  5  20) =#020514     almost black
    background = gf.ocean.opts(facecolor='#01013f') * gf.land.opts(facecolor='#958258') 


t1 = time.time()
plotpoly_hv(xlat,xlon,var,pngname,clim=clim,title=varname,proj=proj,colormap=my_cmap_alpha,width=wres,height=hres,background=background)
t2= time.time()
print(f"hv polycollection:  {t2-t1:.2f}s")


