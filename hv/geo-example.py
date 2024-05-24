import holoviews as hv
import geoviews as gv
import geoviews.feature as gf
from cartopy import crs
#
# modified example from:
# https://geoviews.org/user_guide/Projections.html#
#
gv.extension('matplotlib')

#projplot = crs.Orthographic(central_longitude=-90, central_latitude=30)
projplot = crs.Robinson()
#projplot = crs.PlateCarree())


nyc, beijing = (-74.0, 40.7, 'NYC'), (116.4, 39.9, 'Beijing')
london = (14471.53, 6712008., 'London')

cities_lonlat   = gv.Points([nyc, beijing], vdims='City')
cities_lonlat.opts(color='red')
cities_mercator = gv.Points([london], crs=crs.GOOGLE_MERCATOR, vdims='City')
# convert mercator to lonlat:
cities_lonlat2 = gv.operation.project(cities_mercator, projection=crs.PlateCarree())
cities_lonlat2.opts(color='black')


r = gv.Overlay([gf.coastline,gf.ocean,gf.land]) * cities_lonlat * cities_lonlat2
#r.opts(  gv.opts.Points(color='black'))
r.opts(xlim=(-180.,180))
r.opts(ylim=(-90.,90))
r.opts(data_aspect=1)
r.opts(projection=projplot, global_extent=True)
r.opts(fig_inches=8)

rplot=gv.render(r)
rplot.savefig("geo-example.png", bbox_inches='tight',dpi=600)
rplot.show()
input("pausing to disply plot. hit any key")



