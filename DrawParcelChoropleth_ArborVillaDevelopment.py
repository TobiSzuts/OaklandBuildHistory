# -*- coding: utf-8 -*-
"""
Read a shapefile containing parcel coordinates and housing information, 
and create a choropleth map coded by year each structure was built.

Colorbar definitions and general outline are from 
http://sensitivecities.com/so-youd-like-to-make-a-map-using-python-EN.html#.V2B4hK4tWbj

Focus on Arborvilla

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
from shapely.geometry import Polygon, MultiPolygon, Point
from shapely.prepared import prep
from descartes import PolygonPatch
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
    cmap.set_under = 'magenta'
    cmap.set_over = 'green'
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


# shapefile
baseFile = 'data/Oakland_parcels_queried/Oakland_parcels_queried'


center = geopy.Point(37.8012314, -122.2423577)        # (lat, long), Arbor Villa
radius = 0.4                           # in km
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

# Create projection view as a polygon to filter shapes
window = [ll, (ll[0], ur[1]), ur, (ur[0], ll[1]), ll]
window = list(zip( *m(*list(zip(*window))) ))
window_map = pd.DataFrame({'poly': [Polygon(window)]})
window_polygon = prep(MultiPolygon(list(window_map['poly'].values)))

# Remove any shapes that are outside the map window
df_map = df_map[ [window_polygon.intersects(i) for i in df_map.poly] ]

# draw tract patches from polygons
df_map['patches'] = df_map['poly'].map(lambda x: PolygonPatch(
    x,
    ec='#787878', lw=.25, alpha=.9,
    zorder=4))

# create colormap based on year built
cmap_range = (1924.5, 1951.5)
#cmap_range = (cmap_range[0]-(cmap_range[1]-cmap_range[0])/(ncolors-1), cmap_range[1])
ncolors = 9
yearBuilt_bins = np.linspace(min(cmap_range), max(cmap_range), ncolors+1)
cmap = matplotlib.cm.coolwarm
# set_bad doesn't work with BoundaryNorm
cmap.set_bad(color='white')       # if yearBuilt is nan
cmap.set_under(color=cmap(0))       # if yearBuilt less than min(range)
cmap.set_over(color='white')       # if yearBuilt is more than max(range)
norm = matplotlib.colors.BoundaryNorm(yearBuilt_bins, ncolors)
df_map['color'] = norm(df_map.yearBuilt) 
# normalization changes nans to -1, same as underflow
df_map['color'][df_map.yearBuilt<1700] = -2
#df_map

plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, axisbg='w', frame_on=False)

# plot parcels by adding the PatchCollection to the axes instance
pc = PatchCollection(df_map['patches'].values, match_original=True)
#pc.set_facecolor(cmap(norm(df_map.yearBuilt)/(ncolors-1)))
pc.set_facecolor(cmap(np.ma.masked_less_equal(list(df_map.color), -2)/(ncolors-1)))
ax.add_collection(pc)

yearBuilt_labels = ['%.0f-%.0f' % (np.ceil(yearBuilt_bins[i]), np.floor(yearBuilt_bins[i+1]))
                        for i in range(ncolors)]
#yearBuilt_labels.append('>%.0f' % yearBuilt_bins[-1])
yearBuilt_labels[0] = '<%.0f' % yearBuilt_bins[1]

#yearBuilt_bins = np.insert(yearBuilt_bins, 0, 1700)
#yearBuilt_bins = np.append(yearBuilt_bins, 2016)
cb = colorbar_index(ncolors=ncolors, cmap=cmap, shrink=0.5, 
                    labels=yearBuilt_labels)
#                    boundaries = [i for i in range(ncolors+1)])
#                    , extend='both')
#                    boundaries = [1700] + list(yearBuilt_bins) + [2016])
cb.ax.tick_params(labelsize=6)

# copyright and source data info
smallprint = ax.text(
    1.3, 0.15,
    'Tobi Szuts 2016 \nYear built data from Zillow \nParcel data gathered by Michal Migurski \nin 2011 and downloaded from \nhttp://codeforoakland.org/data-sets/#oakland1',
    ha='right', va='bottom',
    size=4,
    color='#555555',
    transform=ax.transAxes)

for (k, v) in neighborhoods.items() :
    pos = m(v[1], v[0])
    if window_polygon.contains(Point(pos)) is True :
        print(k,v)
        ax.text(pos[0], pos[1],
        k,
        ha='center', va='center',
        size=6,
        color='black')
    
#pos_fidicual = m(-122.2436564, 37.813227)    
#dev = m.scatter(
#    pos_fidicual[0], pos_fidicual[1],
#    5, marker='o', lw=.25,
#    facecolor='#33ccff', edgecolor='w',
#    alpha=0.9, antialiased=True,
#    label='Fiducial', zorder=3)
    
# Draw a map scale
m.drawmapscale(
    coords[0] + w*0.8, coords[1] + h * 0.05,
    coords[0], coords[1],
    1000*radius/2,   # length, in m
    barstyle='fancy', labelstyle='simple',
    units = 'm',
#    format='%.2f',
    fillcolor1='w', fillcolor2='#555555',
    fontcolor='#555555',
    zorder=4)
plt.title("Arbor Villa development, 1926-1935")
plt.tight_layout()
# this will set the image width to 722px at 100dpi
fig.set_size_inches(7.22, 5.25)  # use for larger size
#fig.set_size_inches(3.5, 2.5)
plt.savefig('data/Oakland_temp.png', dpi=300, alpha=True)
plt.show()