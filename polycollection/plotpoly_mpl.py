#
# Matplotlib polycolleciton interface
# plot data on SCRIP polygons, with cartopy
#
#
import numpy as np
from netCDF4 import Dataset
import cartopy 
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy.interpolate import RegularGridInterpolator,interpn

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



def plotpoly(xlat,xlon,data,clat,clon,outname=None, title='',
             proj=ccrs.PlateCarree(), dpi=1200,
             xlim=(-180.,180), ylim=(-90.,90.),
             clim=None,colormap=None,mask=1,colormap_mask=None,
             background='cartopy', interp_bg=False, extent=None
):

    
    # if mask present, remove masked cells
    if not np.isscalar(mask):
        data=data[mask]
        xlon = xlon[mask,:]
        xlat = xlat[mask,:]
        #count = sum(1 for x in mask if x)

    # convert to degrees, if necessary
    if np.max(np.abs(xlat))<1.1*np.pi:
        xlat=xlat*180/np.pi
        xlon=xlon*180/np.pi
        
    mn=float(min(data))
    mx=float(max(data))
    print(f"plotpoly(): {len(data)} cells. data min/max= {mn:.3},{mx:.3} {title}")
    if clim is None:
        clim=(mn,mx)
    if colormap is None:
        if mn*mx < 0: colormap='Spectral'
        else: colormap='plasma'

    # transform into desired coordinate system:
    xpoly  = proj.transform_points(ccrs.PlateCarree(), xlon, xlat)

    # fix and duplicate cut cells
    if "proj=eqc" in proj.srs:
        # duplicate cut polygons on left and right edge of plot
        [xpoly,xpoly_new,mask_new] = shift_anti_meridian_polygons(xpoly)
        corners=np.concatenate((xpoly[:,:,0:2],xpoly_new[:,:,0:2]),axis=0)
        data=np.concatenate((data,data[mask_new]),axis=0)
    if "proj=robin" in proj.srs:
        # remove all cut polygons
        eps=40*1e5
        mask_keep = np.array(np.max(xpoly[:,:,0], axis=1) - np.min(xpoly[:,:,0], axis=1) < eps)
        corners=xpoly[mask_keep,:,0:2]
        data=data[mask_keep]
    if "proj=ortho" in proj.srs:
        #remove non-visible points:
        mask_keep =  np.all(np.isfinite(xpoly),axis=(1,2))
        corners = xpoly[mask_keep,:,0:2]
        data=data[mask_keep]
    if "proj=lcc" in proj.srs:
        #this projection is normally used with an extent
        #polygons which are cut by projection would normally be outside of the extent
        # thus only remove non-visable points, no need to look for cut polygons
        mask_keep = np.all(np.isfinite(xpoly),axis=(1,2))
        corners=xpoly[mask_keep,:,0:2]
        data=data[mask_keep]


    fig=plt.figure(dpi=dpi)                    #create new figure
    ax = plt.axes(projection=proj)      #add axis to the figure
    if extent is None:
        ax.set_global()
    else
        ax.set_extent(extent,crs=ccrs.PlateCarree())  # extent given in lat/lon coords

    # add cartopy background images. needs to be first because it uses facecolor
    background2=None
    if background=='cartopy':
        print("output internal cartopy background")
        # option zorder=0,1,2... will specfy which layer to draw first
        #ax.set_facecolor(cfeature.COLORS['water'])
        ax.set_facecolor('#01013f')
        #ax.coastlines(resolution='110m')
        #ax.add_feature(cartopy.feature.OCEAN,facecolor='#01013f', edgecolor='none')
        #ax.add_feature(cartopy.feature.LAND,edgecolor='black')
        ax.add_feature(cartopy.feature.LAND, zorder=0,facecolor='#958258', edgecolor='none', alpha=1)
        #ax.add_feature(cartopy.feature.LAKES, zorder=100, alpha=0.5)
        ax.add_feature(cartopy.feature.LAKES, zorder=0, facecolor='#01013f', edgecolor='none',alpha=1)
        plt.savefig(f"{outname}-bg.png",dpi='figure',orientation="portrait",bbox_inches='tight',facecolor='white', transparent=False)
        # go back to white background for low-res gaps from missing polygons
        ax.set_facecolor('white')
    
    elif background=='none':
        pass
    else:
        background2 = mpl.image.imread(background)
    
    print("creating polycollection")
#    p = PolyCollection(corners, array=data,edgecolor='face',linewidths=0,antialiased=False)
#    p = PolyCollection(corners, array=data,edgecolor='face',antialiased=False)
    p = PolyCollection(corners, array=data,edgecolor='none',linewidths=0,antialiased=False)

    p.set_clim(clim)
    p.set_cmap(colormap)
    #fig.colorbar(p,ax=ax)
    #ax.set_title(title)

        
    print("add polycollection...")
    ax.add_collection(p)
    # add contenental outlines in black
    #ax.coastlines(resolution='110m') # options: '110m', '50m', '10m'    
    
    print("output polycollection plot...")
    plt.savefig(f"{outname}.png",dpi='figure',orientation="portrait",bbox_inches='tight')

    
    # plot the alpha mask:
    p.set_cmap(colormap_mask)
    print("output polycollection - alpha mask...")
    plt.savefig(f"{outname}-mask.png",dpi='figure',orientation="portrait",bbox_inches='tight')


    if not (background2 is None):
        if interp_bg:
            print("output polycollection background...")
            # pass in xc,yc into here, mask out above
            # interpolate "background" to "rgb"  with coords  linspace() -> xc,yc
            imglon=np.linspace(-180, 180, background2.shape[1])
            imglat=np.linspace(90, -90, background2.shape[0])
            if "proj=eqc" in proj.srs:
                clon2=clon   # mask_keep not computed
                clat2=clat
            else:
                clon2=clon[mask_keep]
                clat2=clat[mask_keep]
                clon2[clon2>180] = clon2[ clon2>180] - 360.
            rgb = interpn( (imglat,imglon), background2,  (clat2,clon2))
            c=rgb[:,0] ;  rgb[c<0,0]=0 ; rgb[c>1,0]=1
            c=rgb[:,1] ;  rgb[c<0,1]=0 ; rgb[c>1,1]=1
            c=rgb[:,2] ;  rgb[c<0,2]=0 ; rgb[c>1,2]=1
            p.set(array=None,facecolors=rgb)
        else:
            print("output imshow background...")
            ax.imshow(background2, zorder=100, extent=[-180, 180, -90, 90], origin='upper', transform=ccrs.PlateCarree())
        plt.savefig(f"{outname}-bg.png",dpi='figure',orientation="portrait",bbox_inches='tight',facecolor='white', transparent=False)


    plt.clf()  # clear figure
    plt.close(fig) # close figure, release memory
