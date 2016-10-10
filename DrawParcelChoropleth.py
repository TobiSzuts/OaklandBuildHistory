# -*- coding: utf-8 -*-
"""
Read a shapefile containing parcel coordinates and housing information, 
and create a choropleth map coded by year each structure was built.

Colorbar definitions and general outline are from 
http://sensitivecities.com/so-youd-like-to-make-a-map-using-python-EN.html#.V2B4hK4tWbj

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
from OaklandNeighborhoodLabels import NeighborhoodLabels

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

# compute map boundscenter = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
center = geopy.Point(37.8058428, -122.2399758)        # (lat, long), Armenian Church
radius = 1.5                           # in km
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
    ec='#787878', lw=.0, alpha=.9,
    zorder=4))
    # lw=0.25 for radius<4 km

# create colormap based on year built
#cmap_range = (1915.5, 1960.5)
#cmap_range = (1885.5, 1930.5)
cmap_range = (1880.5, 2015.5)
ncolors = 9
yearBuilt_bins = np.linspace(min(cmap_range), max(cmap_range), ncolors+1)
cmap = matplotlib.cm.coolwarm
cmap.set_bad(color='white')       # if yearBuilt is nan
cmap.set_under(color=cmap(0))       # if yearBuilt less than min(range)
cmap.set_over(color='red')       # if yearBuilt is more than max(range)
norm = matplotlib.colors.BoundaryNorm(yearBuilt_bins, ncolors)
df_map['color'] = norm(df_map.yearBuilt) 
# normalization changes nans to -1, same as underflow, so remap those
df_map.loc[df_map.yearBuilt<1700, 'color'] = -2

plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, axisbg='w', frame_on=False)

# plot parcels by adding the PatchCollection to the axes instance
pc = PatchCollection(df_map['patches'].values, match_original=True)
pc.set_facecolor(cmap(np.ma.masked_less_equal(list(df_map.color), -2)/(ncolors-1)))
ax.add_collection(pc)

yearBuilt_labels = ['%.0f-%.0f' % (np.ceil(yearBuilt_bins[i]), np.floor(yearBuilt_bins[i+1]))
                        for i in range(ncolors)]
yearBuilt_labels[0] = '<%.0f' % yearBuilt_bins[1]

cb = colorbar_index(ncolors=ncolors, cmap=cmap, shrink=0.5, 
                    labels=yearBuilt_labels)
cb.ax.tick_params(labelsize=6)

# copyright and source data info
smallprint = ax.text(
    1.3, 0.15,
    'Tobi Szuts 2016 \nYear built data from Zillow \nParcel data gathered by Michal Migurski \nin 2011 and downloaded from \nhttp://codeforoakland.org/data-sets/#oakland1',
    ha='right', va='bottom',
    size=4,
    color='#555555',
    transform=ax.transAxes)

#neighborhoods = NeighborhoodLabels()
#for (k, v) in neighborhoods.items() :
#    pos = m(v[1], v[0])
#    if window_polygon.contains(Point(pos)) is True :
#        ax.text(pos[0], pos[1],
#        k,
#        ha='center', va='center',
#        size=6,
#        color='black')

# code to plot a point, used for verifying a new coordinate   
#pos_fidicual = m(-122.2436564, 37.813227)    
#dev = m.scatter(
#    pos_fidicual[0], pos_fidicual[1],
#    5, marker='o', lw=.25,
#    facecolor='#33ccff', edgecolor='w',
#    alpha=0.9, antialiased=True,
#    label='Fiducial', zorder=3)
    
# Draw a map scale
if radius < 1 :
    scaleLen = round(radius*4)/8*1000
else :
    scaleLen = round(radius/2)*1000

m.drawmapscale(
    coords[0] + w*0.8, coords[1] + h * 0.05,
    coords[0], coords[1],
    int(scaleLen),   # length, in m
    barstyle='fancy', labelstyle='simple',
    units = 'm',
    fillcolor1='w', fillcolor2='#555555',
    fontcolor='#555555',
    zorder=4)
plt.title("Oakland housing development, 1920-1960")
plt.tight_layout()
# this will set the image width to 722px at 100dpi
fig.set_size_inches(7.22, 5.25)  # use for larger size
plt.savefig('data/Oakland_temp.png', dpi=300, alpha=True)
plt.show()