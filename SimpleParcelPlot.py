# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 16:09:40 2016

@author: tobis

Script to create a map of immediate neighborhood, no color.
"""


# Load Libraries
#from lxml import etree
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#import matplotlib.cm as cm
#from matplotlib.colors import Normalize
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Point, Polygon, MultiPoint, MultiPolygon
from shapely.prepared import prep
#from pysal.esda.mapclassify import Natural_Breaks as nb
from descartes import PolygonPatch
import fiona
from itertools import chain
import geopy.distance
distance = geopy.distance.vincenty    

# Convenience functions for working with colour ramps and bars
def colorbar_index(ncolors, cmap, labels=None, **kwargs):
    """
    This is a convenience function to stop you making off-by-one errors
    Takes a standard colour ramp, and discretizes it,
    then draws a colour bar with correctly aligned labels
    """
    cmap = cmap_discretize(cmap, ncolors)
    mappable = matplotlib.cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors+0.5)
    colorbar = plt.colorbar(mappable, **kwargs)
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
    colorbar.set_ticklabels(range(ncolors))
    if labels:
        colorbar.set_ticklabels(labels)
    return colorbar
    
def cmap_discretize(cmap, N):
    """
    Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet. 
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)

    """
    if type(cmap) == str:
        cmap = get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki]) for i in range(N + 1)]
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)


baseFile = 'data/alameda-tracts-2000/alameda-tracts-2000'
baseFile = 'data/Oakland_parcels/parcels'
#baseFile = 'data/Oakland_parcels_subset2/parcels'# radius of 1e-3
baseFile = 'data/Oakland_parcels_subset/parcels'# radius of 1e-2

# compute map boundscenter = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
center = geopy.Point(37.8058428, -122.2399758)        # (lat, long), Armenian Church
radius = 0.3                           # in km
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

## Remove any shapes that are outside the map window
#wards_polygon = prep(MultiPolygon(list(df_map['poly'].values)))
#df_map2 = filter(window_polygon.intersects, df_map['poly'])


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
#pc.set_facecolor(cmap(norm(df_map.zip)));
ax.add_collection(pc)
#ax.add_collection(PatchCollection(window_map['patches'].values, match_original=True))

#zip_labels = [str(x+94602) for x in range(10)]
#cb = colorbar_index(ncolors=10, cmap=cmap, shrink=0.5, labels=zip_labels)
#cb.ax.tick_params(labelsize=6)

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
fig.set_size_inches(7.22, 5.25)  # use for larger size
#fig.set_size_inches(3.5, 2.5)
plt.savefig('data/Oakland_temp.png', dpi=600, alpha=True)
plt.show()

