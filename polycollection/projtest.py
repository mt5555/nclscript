import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
#
#  draw contenents on different projections with different extents
#
#  used to quickly choose correct extents for regions of interest
#

pn=4
if pn==1:
    plon=0
    proj=ccrs.PlateCarree(central_longitude=plon) ; projname=f"latlon{plon}"
    NA_extent=[-115,20,0,70]
if pn==2:
    proj=ccrs.Robinson()   ; projname="robinson0"
    NA_extent=[-115,20,0,70]
if pn==3:
    plat=30.; plon=-40.;
    proj = ccrs.Orthographic(central_latitude=plat, central_longitude=plon)
    NA_extent=[plon-50,plon+50,plat-40,plat+35]
if pn==4:
    plon=-45.;
    proj=ccrs.LambertConformal(central_longitude=plon,standard_parallels=(20, 45))
                               #standard_parallels=(33, 45)  default?
                               #standard_parallels=(37, 65)
                               #standard_parallels=(20, 45)  good?
    NA_extent=[plon-50,plon+40,-5,65]
    
print("proj=",proj.srs)
fig=plt.figure(dpi=150)
ax = plt.axes(projection=proj)
#ax.set_extent(NA_extent,crs=ccrs.PlateCarree())
#ax.set_facecolor('#01013f')
#ax.add_feature(cartopy.feature.LAND, zorder=0,facecolor='#958258', edgecolor='none', alpha=1)
#ax.add_feature(cartopy.feature.LAKES, zorder=0, facecolor='#01013f', edgecolor='none',alpha=1)
ax.coastlines()
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=.5, color='black', alpha=0.5, linestyle='-', draw_labels=True)
#gl = ax.gridlines(draw_labels=True)
print("extent (native coords): ",ax.get_extent())

fig=plt.figure(dpi=150)
ax = plt.axes(projection=proj)
ax.set_extent(NA_extent,crs=ccrs.PlateCarree())
ax.set_facecolor('#01013f')
ax.add_feature(cartopy.feature.LAND, zorder=0,facecolor='#958258', edgecolor='none', alpha=1)
ax.add_feature(cartopy.feature.LAKES, zorder=0, facecolor='#01013f', edgecolor='none',alpha=1)
#ax.coastlines()
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=.5, color='white', alpha=0.5, linestyle='-', draw_labels=False)    
    


plt.show()
