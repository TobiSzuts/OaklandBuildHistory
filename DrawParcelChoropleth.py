# -*- coding: utf-8 -*-
"""
Read a shapefile containing parcel coordinates and housing information, 
and create a choropleth map coded by year each structure was built.

Created on Wed Jun  8 09:45:41 2016
@author: tobis
"""

# Load Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, MultiPolygon
from shapely.prepared import prep
from descartes import PolygonPatch
from itertools import chain
import geopy.distance
distance = geopy.distance.vincenty    

# shapefile
baseFile = 'data/Oakland_parcels_queried/Oakland_parcels_queried'


# compute map boundscenter = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
center = geopy.Point(37.8058428, -122.2399758)        # (lat, long), Armenian Church
radius = 0.25                           # in km
ur = distance(kilometers=radius).destination(center, +45)
ll = distance(kilometers=radius).destination(center, -135)
ur = (ur.longitude, ur.latitude)
ll = (ll.longitude, ll.latitude)

extra = 0.01           # padding
coords = list(chain(ll, ur))
w, h = coords[2] - coords[0], coords[3] - coords[1]

m = Basemap(
    projection='tmerc',
    lon_0=-122.,
    lat_0=37.,
    ellps = 'WGS84',
    llcrnrlon=coords[0] - extra * w,
    llcrnrlat=coords[1] - extra * h,
    urcrnrlon=coords[2] + extra * w,
    urcrnrlat=coords[3] + extra * h,
    lat_ts=0,
    resolution='i',
    suppress_ticks=True)

window = [ll, (ll[0], ur[1]), ur, (ur[0], ll[1]), ll]
# convert into map coordinates, unpacking and then repacking tuples
window = list(zip( *m(*list(zip(*window))) ))
window_map = pd.DataFrame({'poly': [Polygon(window)]})
window_polygon = prep(MultiPolygon(list(window_map['poly'].values)))

m.readshapefile(
    baseFile,
    'oakland',
    color='blue',
    zorder=2)
 
  
# set up a map dataframe
df_map = pd.DataFrame({
    'poly': [Polygon(xy) for xy in m.oakland],
    'id': [obj['OBJECTID'] for obj in m.oakland_info],
    'zip': [obj['ZIP'] for obj in m.oakland_info],
    'yearBuilt': [obj['YEARBUILT'] for obj in m.oakland_info]})

# Remove any shapes that are outside the map window
#df_map = df_map[ [window_polygon.intersects(i) for i in df_map.poly] ]

# #%% Plot everything
# draw tract patches from polygons
df_map['patches'] = df_map['poly'].map(lambda x: PolygonPatch(
    x,
    ec='#787878', lw=.25, alpha=.9,
    zorder=4))


# create colormap based on year built
cmap_range = (1880, 1940)
ncolors = 8
yearBuilt_bins = np.linspace(min(cmap_range), max(cmap_range), ncolors+1)
#norm = matplotlib.colors.Normalize(vmin=min(cmap_range), vmax=max(cmap_range))
cmap = matplotlib.cm.coolwarm
cmap.set_bad(color='white')       # if yearBuilt is nan
norm = matplotlib.colors.BoundaryNorm(yearBuilt_bins, ncolors)
#setColor = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, axisbg='w', frame_on=False)

# plot boroughs by adding the PatchCollection to the axes instance
pc = PatchCollection(df_map['patches'].values, match_original=True)
pc.set_facecolor(cmap(norm(df_map.yearBuilt)/ncolors));
ax.add_collection(pc)

yearBuilt_labels = ['%.0f-%.0f' % (yearBuilt_bins[i], yearBuilt_bins[i+1])
                        for i in range(ncolors)]
#yearBuilt_labels[0] = '<%.0f' % yearBuilt_bins[1]
yearBuilt_labels.append('>%.0f' % yearBuilt_bins[-1])

cb = colorbar_index(ncolors=ncolors+1, cmap=cmap, shrink=0.5, labels=yearBuilt_labels)
cb.ax.tick_params(labelsize=6)

# Draw a map scale
m.drawmapscale(
    coords[0] + w * 0.5, coords[1] + h * 0.1,
    coords[0], coords[1],
    radius/2*1000,   # length
    barstyle='fancy', labelstyle='simple',
    units = 'm',
#    format='%.2f',
    fillcolor1='w', fillcolor2='#555555',
    fontcolor='#555555',
    zorder=4)
plt.title("Parcels in Oakland, labeled by year built")
plt.tight_layout()
# this will set the image width to 722px at 100dpi
fig.set_size_inches(7.22, 5.25)  # use for larger size
#fig.set_size_inches(3.5, 2.5)
plt.savefig('data/Oakland_temp.png', dpi=300, alpha=True)
plt.show()