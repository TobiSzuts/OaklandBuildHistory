# -*- coding: utf-8 -*-
"""
Created on Thu May 12 09:56:01 2016

@author: tobis

Further code to develop neighborhood development map

Graphics scratch

"""



#%% Drawing polygons
    
import matplotlib.pyplot as plt
from descartes import PolygonPatch
import shapely.geometry as shp
import geopy
import geopy.distance

data_keys = list(data.keys())       # make into a list to extract only the first couple
data_keys.sort()

fig = plt.figure(1, figsize=(5,5), dpi=100)
ax = fig.add_subplot(111)

#for k in data_keys :
#    e = shp.shape(data[k]['geometry'])
#    f = e.envelope         # make into a simple envelope


e = shp.asShape(data[data_keys[0]]['geometry'])
f_patch = PolygonPatch(e)
ax.add_patch(f_patch)

#e = shp.shape(data[k]['geometry'])
#e = e.envelope

eps = 0.9999
xrange = [min([x for x in e if x<0])*eps, max([x for x in e if x<0])/eps ]
yrange = [min([x for x in e if x>0])*eps, max([x for x in e if x>0])/eps ]

xrange = [-122.22803110589676, -122.2518177180526]
yrange = [37.8017253409998, 37.8098831650553]

ax.set_xlim(*xrange)
ax.set_ylim(*yrange)
#plt.xlim(min(long_disp)*1, max(long_disp))
plt.show()

#%% from http://basemaptutorial.readthedocs.io/en/latest/shapefile.html
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

map = Basemap(llcrnrlon=-0.5,llcrnrlat=39.8,urcrnrlon=4.,urcrnrlat=43.,
             resolution='i', projection='tmerc', lat_0 = 39.5, lon_0 = 2)

map.drawmapboundary(fill_color='aqua')
map.fillcontinents(color='#ddaa66',lake_color='aqua')
map.drawcoastlines()

#map.readshapefile('../sample_files/comarques', 'comarques')

plt.show()

#%% example from https://gist.github.com/urschrei/6436526 to plot random points on 
# a map of London
""" 
required packages:
numpy
matplotlib
basemap: http://matplotlib.org/basemap/users/installing.html
shapely: https://pypi.python.org/pypi/Shapely
descartes: https://pypi.python.org/pypi/descartes
random
numpy and random are only required to generate random points for this example
"""

from random import randint, shuffle
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Point, MultiPoint, MultiPolygon
from descartes import PolygonPatch


# lower left minx miny , upper right maxx maxy
bounds = [-6.108398, 49.61071, 1.669922, 58.972667]
minx, miny, maxx, maxy = bounds
w, h = maxx - minx, maxy - miny

# generate random points within the bounds
lon = np.linspace(minx, maxx).tolist()
lat = np.linspace(miny, maxy).tolist()
random.shuffle(lon)
random.shuffle(lat)

# create a new matplotlib figure and axes instance
fig = plt.figure()
ax = fig.add_subplot(111)
# add a basemap and a small additional extent
m = Basemap(
    projection='merc',
    ellps = 'WGS84',
    llcrnrlon=minx - 0.2 * w,
    llcrnrlat=miny - 0.2 * h,
    urcrnrlon=maxx + 0.2 * w,
    urcrnrlat=maxy + 0.2 * h,
    lat_ts=0,
    resolution='h')
m.drawcoastlines(linewidth=0.3)
m.drawmapboundary()
# a shapefile can be added like so if needed
# m.readshapefile('london_shp', 'london', color='#555555')

# set axes limits to basemap's coordinate reference system 
min_x, min_y = m(minx, miny)
max_x, max_y = m(maxx, maxy)
corr_w, corr_h = max_x - min_x, max_y - min_y
ax.set_xlim(min_x - 0.2 * corr_w, max_x + 0.2 * corr_w)
ax.set_ylim(min_y - 0.2 * corr_h, max_y + 0.2 * corr_h)
# square up axes and basemap
ax.set_aspect(1)
# buffer units are translated to metres by Basemap
# we're randomly varying between 7.5k and 15k metres
patches = [PolygonPatch(Point(m(lon, lat)).buffer(1.0 * randint(7500, 15000)),
                        fc='#cc00cc', ec='#555555', alpha=0.5, zorder=4)
           for lon, lat in zip(lon, lat)]

ax.add_collection(PatchCollection(patches, match_original=True))
#plt.savefig('data/uk.png', dpi=300)
plt.show()


#%% 
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from descartes import PolygonPatch
from shapely.geometry.polygon import Polygon

center = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
# map radius of 1 shows as far south as Santa Cruz, 0.1 has some Bay, 0.01 has no water at all
map_radius = 1e-1

map = Basemap(llcrnrlon=center[1]-map_radius,
              llcrnrlat=center[0]-map_radius,
              urcrnrlon=center[1]+map_radius,
              urcrnrlat=center[0]+map_radius, 
              resolution = 'h', projection='aeqd', lat_0=center[0],
              lon_0 = center[1]);
    # 'cyl' = equidistant cylindrical projection
    # 'aeqd' = Azimuthal Equidistant Projection              
              
fig = plt.figure()
ax = fig.add_subplot(111)

# set axis limits to basemap coordinate system
#ax.set_xlim(center[1]-map_radius, center[1]+map_radius)
#ax.set_ylim(center[0]-map_radius, center[0]+map_radius)
#ax.set_aspect(1)        # square up axes

map.drawmapboundary(fill_color='aqua')
#map.fillcontinents(color='#ddaa66',lake_color='aqua')
map.drawcoastlines()

lon = [-122.2399758, -122.2399758, -122.2399758 + 0.01, -122.2399758 + 0.01, -122.2399758]
lat = [37.8058428, 37.8058428+0.01, 37.8058428+0.01, 37.8058428, 37.8058428]

lon = [-122.2399758]
lat = [37.8058428]

#patches = [PolygonPatch(Point(map(lon, lat)).buffer(1.0 * randint(7500, 15000)),
#                        fc='#cc00cc', ec='#555555', alpha=0.5, zorder=2)
#           for lon, lat in zip(lon, lat)]
patches = [PolygonPatch(Point(map(lon, lat)),
                        fc='#cc00cc', ec='#555555', alpha=0.5, zorder=6)
           for lon, lat in zip(lon, lat)]
ax.add_collection(PatchCollection(patches, match_original=True))
              


#map.readshapefile('../sample_files/comarques', 'comarques')

plt.show()    

#%% New tack, working off of web site
# http://sensitivecities.com/so-youd-like-to-make-a-map-using-python-EN.html#.V0Xi6VrzQ2w
# DOESN'T WORK!!

import fiona
from itertools import chain
from mpl_toolkits.basemap import Basemap

#shp = fiona.open('data/AlamedaCounty_StreetCenterlines/streets.shp')
shp = fiona.open('data/alameda-openstreetmap/alameda-openstreetmap-polygons.shp')
bds = shp.bounds
shp.close()

coords = list(bds)
w, h = coords[2] - coords[0], coords[3] - coords[1]

center = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
extra = 0.1

m = Basemap(llcrnrlon=coords[0] - extra*w,
              llcrnrlat=coords[1] - extra + 0.01*h,
              urcrnrlon=coords[2] + extra*2,
              urcrnrlat=coords[3] + extra + 0.01*h, 
              resolution = 'i',
              projection='aeqd',
              lat_0=coords[1]+w/2, 
#              suppress_ticks=True,
              lon_0=coords[0]+h/2);
    # 'cyl' = equidistant cylindrical projection
    # 'aeqd' = Azimuthal Equidistant Projection  

m.readshapefile(
    'data/alameda-openstreetmap/alameda-openstreetmap-polygons',
    'oakland',
    color='blue',
    zorder=2)

import pandas as pd
from shapely.geometry import Point, Polygon, MultiPoint, MultiPolygon
from descartes import PolygonPatch
from shapely.prepared import prep

# set up map dataframe
df_map = pd.DataFrame({'poly': [Polygon(x) for x in m.oakland], 
                       'street_name': [s['NAME'] for s in m.oakland_info]})
wards_polygon = prep(MultiPolygon(list(df_map['poly'].values)))

                      
df_map['patches'] = df_map['poly'].map(lambda x: descartes.PolygonPatch(
                                        x, fc='#555555',
                                        ec='#787878', lw=0.25,
                                        alpha=0.9,
                                        zorder=4))

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection

plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, axisbg='w', frame_on=False)

# add shape collection
ax.add_collection(PatchCollection(df_map['patches'].values, match_original=True))

plt.title("Oakland streets")
plt.tight_layout()
#fig.set_size_inches(7)
plt.savefig('TestPlot.png', dpi=100, alpha=True)
plt.show()

#%% Another attempt - works, Modified for Oakland

#%% Load Libraries
from lxml import etree
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Point, Polygon, MultiPoint, MultiPolygon
from shapely.prepared import prep
from pysal.esda.mapclassify import Natural_Breaks as nb
from descartes import PolygonPatch
import fiona
from itertools import chain


baseFile = 'data/Oakland_parcels/parcels'

shp = fiona.open(baseFile + '.shp')
bds = shp.bounds
shp.close()
extra = 0.01
ll = (bds[0], bds[1])
ur = (bds[2], bds[3])
coords = list(chain(ll, ur))
w, h = coords[2] - coords[0], coords[3] - coords[1]

m = Basemap(
    projection='tmerc',
    lon_0=-122.,
    lat_0=37.,
    ellps = 'WGS84',
    llcrnrlon=coords[0] - extra * w,
    llcrnrlat=coords[1] - extra + 0.01 * h,
    urcrnrlon=coords[2] + extra * w,
    urcrnrlat=coords[3] + extra + 0.01 * h,
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
    'ward_name': [ward['OBJECTID'] for ward in m.oakland_info]})
#df_map['area_m'] = df_map['poly'].map(lambda x: x.area)
#df_map['area_km'] = df_map['area_m'] / 100000

## Create Point objects in map coordinates from dataframe lon and lat values
wards_polygon = prep(MultiPolygon(list(df_map['poly'].values)))

# #%% Plot everything
# draw tract patches from polygons
df_map['patches'] = df_map['poly'].map(lambda x: PolygonPatch(
    x,
    fc='#555555',
    ec='#787878', lw=.25, alpha=.9,
    zorder=4))

plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, axisbg='w', frame_on=False)

# plot boroughs by adding the PatchCollection to the axes instance
ax.add_collection(PatchCollection(df_map['patches'].values, match_original=True))
# copyright and source data info
#smallprint = ax.text(
#    1.03, 0,
#    'Total points: %s\nContains Ordnance Survey data\n$\copyright$ Crown copyright and database right 2013\nPlaque data from http://openplaques.org' % len(ldn_points),
#    ha='right', va='bottom',
#    size=4,
#    color='#555555',
#    transform=ax.transAxes)

# Draw a map scale
m.drawmapscale(
    coords[0] + 0.08, coords[1] + 0.015,
    coords[0], coords[1],
    10.,
    barstyle='fancy', labelstyle='simple',
    fillcolor1='w', fillcolor2='#555555',
    fontcolor='#555555',
    zorder=5)
plt.title("Tax Parcels in Oakland")
plt.tight_layout()
# this will set the image width to 722px at 100dpi
fig.set_size_inches(7.22, 5.25)
plt.savefig('data/Oakland_temp.png', dpi=300, alpha=True)
plt.show()


#%% Read shape file, filter, and save a smaller size

import fiona
baseFile = 'data/Oakland_parcels/parcels'

center = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
radius = 1e-2                           # should be 0.001
ll = (center[1]-radius*1, center[0]-radius)
ur = (center[1]+radius*1, center[0]+radius)
ROI = (*ll, *ur)

counter = 0
# Register format drivers with a context manager
with fiona.drivers():

    # Open a file for reading. We'll call this the "source."
    with fiona.open(baseFile + '.shp') as source:

        # The file we'll write to, the "sink", must be initialized
        # with a coordinate system, a format driver name, and
        # a record schema.  We can get initial values from the open
        # collection's ``meta`` property and then modify them as
        # desired.

        meta = source.meta

        # Open an output file, using the same format driver and
        # coordinate reference system as the source. The ``meta``
        # mapping fills in the keyword parameters of fiona.open().

        with fiona.open('data/Oakland_parcels_subset_2/parcels.shp', 'w', **meta) as sink:

            # Process only the records intersecting a box.
            for f in source.filter(bbox=ROI) :
                sink.write(f)
                counter += 1

    # The sink's contents are flushed to disk and the file is
    # closed when its ``with`` block ends. This effectively
    # executes ``sink.flush(); sink.close()``.

# At the end of the ``with fiona.drivers()`` block, context
# manager exits and all drivers are de-registered.

#%% Now modify to show a smaller extent of Oakland, by filtering for all objects inside map frame
#%% Another attempt - works, Modified for Oakland

# Load Libraries
from lxml import etree
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Point, Polygon, MultiPoint, MultiPolygon
from shapely.prepared import prep
from pysal.esda.mapclassify import Natural_Breaks as nb
from descartes import PolygonPatch
import fiona
from itertools import chain
import geopy.distance
distance = geopy.distance.vincenty    



baseFile = 'data/alameda-tracts-2000/alameda-tracts-2000'
baseFile = 'data/Oakland_parcels/parcels'
baseFile = 'data/Oakland_parcels_subset2/parcels'# radius of 1e-3
baseFile = 'data/Oakland_parcels_subset/parcels'# radius of 1e-2

# compute map bounds
center = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
radius = 1e-3                           # should be 0.001

#---
origin = geopy.Point(lat1, lon1)
destination = VincentyDistance(kilometers=d).destination(origin, b)
lat2, lon2 = destination.latitude, destination.longitude
#---

ll = (center[1]-radius*10, center[0]-radius)
ur = (center[1]+radius*10, center[0]+radius)


#shp = fiona.open(baseFile + '.shp')
#bds = shp.bounds
#shp.close()
extra = 0.001
#ll = (bds[0], bds[1])
#ur = (bds[2], bds[3])
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
    'ward_name': [obj['OBJECTID'] for obj in m.oakland_info]})
#df_map['area_m'] = df_map['poly'].map(lambda x: x.area)
#df_map['area_km'] = df_map['area_m'] / 100000

## Remove any shapes that are outside the map window


#wards_polygon = prep(MultiPolygon(list(df_map['poly'].values)))
#
#df_map2 = filter(window_polygon.intersects, df_map['poly'])


# #%% Plot everything
# draw tract patches from polygons
df_map['patches'] = df_map['poly'].map(lambda x: PolygonPatch(
    x,
    fc='#555555',
    ec='#787878', lw=.25, alpha=.9,
    zorder=4))
#window_map['patches'] = window_map['poly'].map(lambda x: PolygonPatch(
#    x,
#    fc='blue',
#    ec='blue', lw=.25, alpha=.9,
#    zorder=5))


plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, axisbg='w', frame_on=False)

# plot boroughs by adding the PatchCollection to the axes instance
ax.add_collection(PatchCollection(df_map['patches'].values, match_original=True))
#ax.add_collection(PatchCollection(window_map['patches'].values, match_original=True))

#m.scatter(m(center)


#m.scatter(
#    [geom.x for geom in ldn_points],
#    [geom.y for geom in ldn_points],
#    5, marker='o', lw=.25,
#    facecolor='#33ccff', edgecolor='w',
#    alpha=0.9, antialiased=True,
#    label='Blue Plaque Locations', zorder=3)

# copyright and source data info
#smallprint = ax.text(
#    1.03, 0,
#    'Total points: %s\nContains Ordnance Survey data\n$\copyright$ Crown copyright and database right 2013\nPlaque data from http://openplaques.org' % len(ldn_points),
#    ha='right', va='bottom',
#    size=4,
#    color='#555555',
#    transform=ax.transAxes)

# Draw a map scale
# find width of map in km


m.drawmapscale(
    coords[0] + w * 0.5, coords[1] + h * 0.1,
    coords[0], coords[1],
    0.5,   # length
    barstyle='fancy', labelstyle='simple',
    fillcolor1='w', fillcolor2='#555555',
    fontcolor='#555555',
    zorder=5)
plt.title("Parcels in Oakland")
plt.tight_layout()
# this will set the image width to 722px at 100dpi
fig.set_size_inches(7.22, 5.25)  # use for larger size
#fig.set_size_inches(3.5, 2.5)
plt.savefig('data/Oakland_temp.png', dpi=300, alpha=True)
plt.show()

#%% Latest attempt, 5/31, improved
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
radius = 0.5                           # in km
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
    fc='#555555',
    ec='#787878', lw=.25, alpha=.9,
    zorder=4))

# create colormap based on zip codes, from 94602-94611
norm = matplotlib.colors.Normalize(vmin=94602, vmax=94611)
cmap = matplotlib.cm.Set3
#setColor = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, axisbg='w', frame_on=False)

# plot boroughs by adding the PatchCollection to the axes instance
pc = PatchCollection(df_map['patches'].values, match_original=True)
    # 'match_original' = keep color of original patch
pc.set_facecolor(cmap(norm(df_map.zip)));
ax.add_collection(pc)
#ax.add_collection(PatchCollection(window_map['patches'].values, match_original=True))

zip_labels = [str(x+94602) for x in range(10)]
cb = colorbar_index(ncolors=10, cmap=cmap, shrink=0.5, labels=zip_labels)
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
plt.title("Parcels in Oakland, labeled by zipcode")
plt.tight_layout()
# this will set the image width to 722px at 100dpi
fig.set_size_inches(7.22, 5.25)  # use for larger size
#fig.set_size_inches(3.5, 2.5)
plt.savefig('data/Oakland_temp.png', dpi=300, alpha=True)
plt.show()
