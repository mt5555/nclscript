import numpy as np
#
import time
import sys
import os
from netCDF4 import Dataset
from pathlib import Path
import pickle


import matplotlib as mpl
import geoviews as gv
import geoviews.feature as gf
from plotpoly_mpl import plotpoly
from plotpoly_hv import plotpoly as plotpoly_hv
import cartopy.crs as ccrs
from PIL import Image   # needed to load JPG background image


#
# notes on background image processing:
# see README for details and motivation
#
# interp_bg = False:  use "imshow" to draw image 
# interp_bg = True:
#    MPL: intepolate background image to data polygons
#    HV:  project background image before plotting (builds new triagnularizaiton, SLOW)
# 
# data <= NE120    interp_bg=False (imshow)     5400x2700 image
# data < NE1024    interp_bg=True               5400x2700 image
# data >= NE1024   interp_bg=True               21600x10800 image
#



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


# location of backround images and HV pickle files
bg_path = os.path.expanduser('~/scratch1/viz/bg')

#name="/home/ac.mtaylor/scratch1/mapping/grids/TEMPEST_ne30pg2.scrip.nc"
#dname="/home/ac.mtaylor/scratch1/viz/1995-decadal/ne30-all.nc"
#image_path=f"{bg_path}/world.topo.200408.3x5400x2700.png"
#interp_bg=False
#idx=-1  # last frame


name="/home/ac.mtaylor/scratch1/mapping/grids/TEMPEST_ne1024pg2.scrip.nc"
dname="/home/ac.mtaylor/scratch1/viz/1995-testb2/data/output.scream.decadal.1hourlyINST_ne1024pg2.INSTANT.nhours_x1.1995-08-31-03600.nc"
image_path=f"{bg_path}/world.topo.200408.3x21600x10800.png"
interp_bg=True
idx=2  # 3rd frame, time=334.125
idx=17 # 334.79     time=334.75  

if len(sys.argv)==4:
    name=sys.argv[1]
    dname=sys.argv[2]
    bgopt=int(sys.argv[3])
    if bgopt == 1:   # suggested for resolution up to NE120(27km)
                     # imshow low-res backgroudn
        interp_bg=False
        image_path=f"{bg_path}/world.topo.200408.3x5400x2700.png"
    elif bgopt == 2:  # suggested for resolutions  120 < NE < 1024
        interp_bg=True
        image_path=f"{bg_path}/world.topo.200408.3x5400x2700.png"
    else:            # best results, but slow
        interp_bg=True
        image_path=f"{bg_path}/world.topo.200408.3x21600x10800.png"
    idx=None
elif len(sys.argv)==1:
    # sent vars by hand above
    pass    
else:
    print('drivery.py requires 3 arguments:')
    print('driver.py  scripfile  datafile  resolution_in_km')
    sys.exit();




print(f"input: {dname}")
print(f"background image: {image_path} interp={interp_bg}")
file1 = Dataset(name,"r")
ncols = file1.dimensions["grid_size"].size
clat = file1.variables["grid_center_lat"][:]
clon = file1.variables["grid_center_lon"][:]
xlat = file1.variables["grid_corner_lat"][:,:]
xlon = file1.variables["grid_corner_lon"][:,:]
#area = file1.variables["grid_area"][:]
#Rearth_km = 6378.1                # radius of earth, in km
#var=Rearth_km*np.sqrt(area)

datafile=Dataset(dname,"r")
#varname="Preciptable Water"   ;  varnamef="VapWaterPath"  ; pngname='tmq'
varname="LW"   ;  varnamef="LW_flux_up_at_model_top"      ; pngname='lw'
#varname="SW"   ;  varnamef="SW_flux_up_at_model_top"       ; pngname='sw'
#varname="Vapor (2m)"   ;  varnamef="qv_2m"                ; pngname='qv2m'
#varname="Precip"   ;  varnamef="precip_total_surf_mass_flux"  ; pngname='prec'

dtime_all = datafile.variables["time"][:]  # times
print("times min,max=",np.min(dtime_all),np.max(dtime_all))

pn=3
extent=None  # use global, unless specified below
background_is_fixed = True
if pn==1:
    plon=0
    proj=ccrs.PlateCarree(central_longitude=plon) ; projname=f"latlon{plon}"
    # if projection exactly matches background image, use imshow, dont interpolate
    if plon==0:
        interp_bg=False
    wres=8000 ; hres=round(wres/2)   
    dpi=1600
if pn==2:
    proj=ccrs.Robinson()   ; projname="robinson0"
    wres=10000 ; hres=round(wres/2)
    dpi=1600    # mpl image: 8.3K x 4.3K     (NE1024: 8k pts on equator)
if pn==3:
    plat=28.; plon=-53.;   # atlantic
    plat=28.; plon=-85     #shifted to show some pacific
    #plon = 180. - np.mod(dtime,1)*360.   # follow the sun
    proj = ccrs.Orthographic(central_latitude=plat, central_longitude=plon)
    #extent=[plon-51.5,plon+51.5,plat-32,plat+32]
    projname=f"ortho_{plat:.0f}_{plon:.0f}"
    wres=2000 ; hres=wres    # ne1024/ortho needs wres>2000 to avoid speckling
    dpi=1000   #             # ne1024  4K pts visable, should for 4K x 4K image
if pn==4:  # good for zoomed in over CONUS
    plon=-45.;
    proj=ccrs.LambertConformal(central_longitude=plon,standard_parallels=(20, 45))
                               #standard_parallels=(33, 45)  # default
    #v1: extent=[plon-50,plon+40,-5,65]
    extent=[plon-50,plon+32,5,50]
    projname=f"lcc-na1"
    wres=500 ; hres=wres    # ne1024/ortho needs wres>2000 to avoid speckling
    dpi=400                 # ne1024  4K pts visable, should for 4K x 4K image

print(proj.srs)





# create 3 colormaps:
#    my_cmap         =  MPL: standard, without transparency
#    my_cmap_alpha   =  HV: colormap with alpha channel
#    my_cmap_mask    =  MPL: colormap that will produce an image with alpha values
#
if varname=="SW":
    clim=(0.,1200.)
    cmap = mpl.pyplot.get_cmap("Greys_r")(np.linspace(0, 1, 256))
    my_cmap = mpl.colors.ListedColormap(cmap)

    cmap=cmap.copy()  #  make a new copy otherwise we change my_cmap
    # 100% white, with transparancy
    cmap[:,0]=1 ; cmap[:,1]=1 ; cmap[:,2]=1

    ramp_x0=50     # start around 50 so land is bright/clear 
    ramp_x=64      # stopping around 90 dims the clouds too much
    ramp_y=1.0
    cmap[0:ramp_x0]=0
    #cmap[0,3]=1.0      # opaque where SW=0 (dark)
                        # beter to use ortho projection and dont show dark regions
    cmap[ramp_x0:ramp_x,3]=np.linspace(0,ramp_y,ramp_x-ramp_x0)
    cmap[ramp_x:256,3]=np.linspace(ramp_y,1,256-ramp_x)
    my_cmap_alpha = mpl.colors.ListedColormap(cmap)
    alpha=cmap[:,3]



if varname=="LW":
    clim=(65.,375.)
    cmap = mpl.pyplot.get_cmap("Greys")(np.linspace(0, 1, 256))
    my_cmap = mpl.colors.ListedColormap(cmap)

    cmap=cmap.copy()  #  make a new copy otherwise we change my_cmap
    # 100% white, with transparancy
    cmap[:,0]=1 ; cmap[:,1]=1 ; cmap[:,2]=1

    ramp_x0=126
    ramp_x =200
                    
    ramp_y=1.0
    cmap[0:ramp_x0]=1
    cmap[ramp_x0:ramp_x,3]=np.linspace(ramp_y,0,ramp_x-ramp_x0)
    cmap[ramp_x:256,3]=0
    my_cmap_alpha = mpl.colors.ListedColormap(cmap)
    alpha=cmap[:,3]



    

# GOOD SETTINGS FOR TMQ:
if varname=="Preciptable Water":
    cmap = mpl.pyplot.get_cmap("RdBu_r")(np.linspace(0.1, .95, 256))
    my_cmap = mpl.colors.ListedColormap(cmap)   

    cmap=cmap.copy()  #  make a new copy otherwise we change my_cmap
    clim=(0.,90.)
    opaque=round(.55*my_cmap.N)

    if extent is not None:  # zoomed in regions, less orange, more transparancy
        clim=(0.,110)
        #clim=(0.,100) also good
        opaque=round(.65*my_cmap.N)
    
    cmap[0:opaque, 3] = np.linspace(0, 1, opaque)
    my_cmap_alpha = mpl.colors.ListedColormap(cmap)
    alpha=cmap[:,3]



# create a colormap with colors = alpha 
cmap=cmap.copy()  #  make a new copy otherwise we change my_cmap
cmap[:,0]=alpha
cmap[:,1]=alpha
cmap[:,2]=alpha
cmap[:,3]=1
my_cmap_mask = mpl.colors.ListedColormap(cmap,'mt3')










if idx is None:
    idx1=0
    idx2=len(dtime_all)
else:
    idx1=idx
    idx2=idx+1


#
#  MPL plot
#

if True:
    background = image_path
    #background = 'cartopy'   #dark blue ocean, brown land

    for idx in range(idx1,idx2):
        dtime=dtime_all[idx]
        t0= time.time()
        print("reading data at time=",dtime,flush=True)
        var = datafile.variables[varnamef][idx,:]  # last snapshot in file
        #print("nan_count =",np.count_nonzero(np.isnan(var)))
        pngname2 = f"{pngname}-{dtime:.2f}"

        t1= time.time()
        print(f"data read time: {t1-t0:.2f}s",flush=True)

        if pn==3 and  background_is_fixed == False:
            plon = 180. - np.mod(dtime,1)*360.   # follow the sun
            proj = ccrs.Orthographic(central_latitude=plat, central_longitude=plon)

        plotpoly(xlat,xlon,var,clat,clon,pngname2,title=varname,proj=proj,colormap=my_cmap,colormap_mask=my_cmap_mask,
                 clim=clim,dpi=dpi,background=background,interp_bg=interp_bg,extent=extent)

        t2= time.time()
        print(f"mpl polycollection: {t2-t1:.2f}s",flush=True)

        if background != 'none':
            bg_out_name=f"{pngname2}-bg.png"
        if background_is_fixed:
            background = 'none'    # skip recomputing background for future plotpoly()

        
        print("running magick composite...",bg_out_name,flush=True)
        if 0==os.system(f"magick  {pngname2}-mask.png -colorspace Gray -flatten {pngname2}-mask2.png"):
            os.system(f"magick composite {pngname2}.png {bg_out_name} {pngname2}-mask2.png {pngname2}-composite.png")
        else:
            print("Error running magick to composit image and background")


        t3= time.time()
        print(f"magick composit time: {t3-t2:.2f}s",flush=True)        

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
    projname=f"{bg_path}/world-hr-{projname}.pickle"
    p=Path(projname)
    if p.exists() and p.is_file():
        print(f"loading background image {p}")
        with open(projname, 'rb') as file:
            background = pickle.load(file)
    else:
        Image.MAX_IMAGE_PIXELS = 233280000
        img_data = np.array(Image.open(image_path)) 
        bounds = (-180, 90, 180, -90)  # Assuming the image covers the whole globe
        background=gv.RGB((np.linspace(-180, 180, img_data.shape[1]),
                           np.linspace(90, -90, img_data.shape[0]),
                           img_data[..., 0], img_data[..., 1], img_data[..., 2]),
                          bounds=bounds,crs=ccrs.PlateCarree())  #.opts(projection=proj)
        if interp_bg:
            print(f"projecting background image {p} SLOW!")
            background = gv.project(background,projection=proj)
            print("saving projected background for future use")
            with open(projname, 'wb') as file:
                pickle.dump(background, file)
        else:
            # for latlon projections, dont need to project:
            # what about rotated latlon???
            background.opts(projection=proj)
if bg==2:
    background = gv.Overlay([gf.coastline,gf.ocean,gf.land])
if bg==3:
    #https://www.color-hex.com/color-palette/1021516  ocean: #004589  land: #72601b 
    # blue from 2012 ANL video:  rgb(48 62 141) = #303E8D   nice blue
    # blue from me:  #01013f                                dark blue
    # blue from world iamge      rgb(2  5  20) =#020514     almost black
    background = gf.ocean.opts(facecolor='#01013f') * gf.land.opts(facecolor='#958258') 


idx=-1  # last image in file
dtime=dtime_all[idx]
print("reading data at time=",dtime)
var = datafile.variables[varnamef][idx,:]  
pngname2 = f"{pngname}-{dtime:.2f}"    
    
t1 = time.time()
plotpoly_hv(xlat,xlon,var,pngname2,clim=clim,title=varname,proj=proj,colormap=my_cmap_alpha,width=wres,height=hres,background=background)
t2= time.time()
print(f"hv polycollection:  {t2-t1:.2f}s")


