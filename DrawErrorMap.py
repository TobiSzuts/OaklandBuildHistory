# -*- coding: utf-8 -*-
"""
Read a shapefile containing parcel coordinates and housing information, 
and create a map showing parcel boundaries and error addresses.

Created on Mon Jun 13 10:30:38 2016
@author: tobis
"""

# Load Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, MultiPolygon, Point
from shapely.prepared import prep
from descartes import PolygonPatch
from itertools import chain
import geopy.distance
distance = geopy.distance.vincenty    

# shapefile
baseFile = 'data/Oakland_parcels_queried/Oakland_parcels_queried'
errorFile = 'data/Oakland_parcels_errors/Oakland_parcels_errors'

# compute map boundscenter = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
#center = geopy.Point(37.8058428, -122.2399758)        # (lat, long), Armenian Church
center = geopy.Point(37.79539889, -122.21547850)        # (lat, long), rough geometric center of entire dataset
#center = geopy.Point(37.759987, -122.221965)        # (lat, long), rough geometric center of error dataset
radius = 8.5                       # in km
#ur = distance(kilometers=radius).destination(center, +45)
#ll = distance(kilometers=radius).destination(center, -135)
ur = distance(kilometers=radius*2**0.5).destination(center, +45)
ll = distance(kilometers=radius*2**0.5).destination(center, -135)
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
# draw tract patches from polygons
df_map['patches'] = df_map['poly'].map(lambda x: PolygonPatch(
    x, fc='none',
    ec='#888888', lw=.25, alpha=.9,
    zorder=4))

# color code the errors
for (yb, p, i) in zip(df_map.yearBuilt, df_map.patches, range(50000)) :
    if np.isnan(yb) :
        p.set_color('green')
    elif yb<1830 :          # because the oldest building in Oakland is from the 1850s
        p.set_color('blue')
    elif yb>2016 :
        p.set_color('magenta')
    
## create colormap based on year built
#cmap_range = (1880, 1940)
#ncolors = 8
#yearBuilt_bins = np.linspace(min(cmap_range), max(cmap_range), ncolors+1)
#cmap = matplotlib.cm.coolwarm
#cmap.set_bad(color='white')       # if yearBuilt is nan
#norm = matplotlib.colors.BoundaryNorm(yearBuilt_bins, ncolors)

m.readshapefile(
    errorFile,
    'error',
    color='red',
    zorder=2)

# set up a map dataframe
df_map_error = pd.DataFrame({
    'poly': [Polygon(xy) for xy in m.error],
    'id': [obj['OBJECTID'] for obj in m.error_info],
    'zip': [obj['ZIP'] for obj in m.error_info],
    'yearBuilt': [obj['YEARBUILT'] for obj in m.error_info]})
# draw tract patches from polygons
df_map_error['patches'] = df_map_error['poly'].map(lambda x: PolygonPatch(
    x, fc='red',
    ec='black', lw=.25, alpha=.9,
    zorder=4))

plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, axisbg='w', frame_on=False)

# plot boroughs by adding the PatchCollection to the axes instance
ax.add_collection(PatchCollection(df_map['patches'].values, match_original=True))
ax.add_collection(PatchCollection(df_map_error['patches'].values, match_original=True))

## draw errors as points
#dev = m.scatter([p.x for p in data_errors_points], [p.y for p in data_errors_points],
#          5, marker='o', lw=.25,
#          facecolor='red', edgecolor='black',
#          alpha=0.9, antialiased=True,
#          label='Error Locations')

ncolors = 2
cmap = matplotlib.colors.ListedColormap(['b', 'r'])
norm = matplotlib.colors.BoundaryNorm([0, 10, 20], ncolors)
error_labels = ['Invalid date', 'Address not \nin Zillow']
cb = colorbar_index(ncolors=ncolors, cmap=cmap, shrink=0.25, 
                    labels=error_labels, aspect=10)

cb.ax.tick_params(labelsize=6)
cb.ax.tick_params(size=0)

if False :
    neighborhoods = NeighborhoodLabels()
    for (k, v) in neighborhoods.items() :
        pos = m(v[1], v[0])
        if window_polygon.contains(Point(pos)) is True :
            ax.text(pos[0], pos[1],
            k,
            ha='center', va='center',
            size=6,
            color='black')

if radius < 2 :
    scaleLen = round(radius*4)/8*1000
else :
    scaleLen = round(radius/2)*1000

# Draw a map scale
m.drawmapscale(
    coords[0] + w * 0.2, coords[1] + h * 0.1,
    coords[0], coords[1],
    int(scaleLen),   # length in meters
    barstyle='fancy', labelstyle='simple',
    units = 'm',
    fillcolor1='w', fillcolor2='#555555',
    fontcolor='#555555',
    zorder=4)
plt.title("Parcels in Oakland with an error")
plt.tight_layout()
# this will set the image width to 722px at 100dpi
fig.set_size_inches(7.22, 5.25)  # use for larger size
plt.savefig('data/Oakland_temp.png', dpi=300, alpha=True)
plt.show()