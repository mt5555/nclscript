import holoviews as hv
import geoviews as gv
import geoviews.feature as gf
import matplotlib.pyplot as plt
from bokeh.plotting import show

import matplotlib
#matplotlib.use('Qt5Agg')
 
from geoviews import opts
from cartopy import crs


#gv.extension('bokeh', 'matplotlib')
gv.extension('matplotlib')


nyc, beijing = (-74.0, 40.7, 'NYC'), (116.4, 39.9, 'Beijing')
london = (14471.53, 6712008., 'London')

cities_lonlat   = gv.Points([nyc, beijing], vdims='City')
cities_mercator = gv.Points([london], crs=crs.GOOGLE_MERCATOR, vdims='City')

#(gv.tile_sources.OSM * cities_lonlat * cities_mercator).opts(
#    opts.Points(global_extent=True, width=500, height=475, size=12, color='black'))


#plot = gv.tile_sources.OSM * cities_lonlat * cities_mercator
#plot.opts(  opts.Points(global_extent=True, width=500, height=475, size=12, color='black'))

#features = gv.Overlay([gf.ocean, gf.land, gf.rivers, gf.lakes, gf.borders, gf.coastline])
#gv.output(features, backend='matplotlib', fig='png', size=300)

r = gv.Overlay([gf.coastline]) * cities_lonlat 
r.opts(xlim=(-180.,180))
r.opts(ylim=(-90.,90))
r.opts(data_aspect=1)
#r.opts(frame_width=100)

rplot=gv.render(r,backend='matplotlib')
#rplot.savefig("geo-example.png", bbox_inches='tight',dpi=300)

rplot.show()
input("pausing to disply plot. hit any key")



