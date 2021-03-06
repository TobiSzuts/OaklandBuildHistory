# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 16:09:40 2016

@author: tobis

Script to create a map of immediate neighborhood, no color, black outlines only.
"""


# Load Libraries
#from lxml import etree
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
from shapely.prepared import prep
from descartes import PolygonPatch
from itertools import chain
import geopy.distance
distance = geopy.distance.vincenty    

baseFile = 'data/Oakland_parcels_subset/parcels'# radius of 1e-2

# compute map boundscenter = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
center = geopy.Point(37.8058428, -122.2399758)        # (lat, long), Armenian Church
radius = 0.4                           # in km
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
    'yearBuilt': [1950 for obj in m.oakland_info]})

# #%% Plot everything
# draw tract patches from polygons
df_map['patches'] = df_map['poly'].map(lambda x: PolygonPatch(
    x,
    fc='#FFFFFF',
    ec='#000000', lw=.25, alpha=.9,
    zorder=4))

plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, axisbg='w', frame_on=False)

# plot boroughs by adding the PatchCollection to the axes instance
pc = PatchCollection(df_map['patches'].values, match_original=True)
    # 'match_original' = keep color of original patch
ax.add_collection(pc)

# Draw a map scale
m.drawmapscale(
    coords[0] + w * 0.8, coords[1] + h * 0.04,
    coords[0], coords[1],
    radius/2*1000,   # length
    barstyle='fancy', labelstyle='simple',
    units = 'm',
#    format='%.2f',
    fillcolor1='w', fillcolor2='#555555',
    fontcolor='#555555',
    zorder=4)
plt.title("Parcels in Oakland")
plt.tight_layout()
# this will set the image width to 722px at 100dpi
fig.set_size_inches(7.5, 7.5)  # use for larger size
#plt.savefig('data/Oakland_ClevelandHeights_outline.png', dpi=600, alpha=True)
plt.savefig('data/Oakland_ClevelandHeights_outline2.pdf', dpi=600, alpha=True)
plt.show()

