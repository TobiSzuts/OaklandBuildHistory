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
map_radius = 1e-0

map = Basemap(llcrnrlon=center[1]-map_radius,
              llcrnrlat=center[0]-map_radius,
              urcrnrlon=center[1]+map_radius,
              urcrnrlat=center[0]+map_radius, 
              resolution = 'f', projection='aeqd', lat_0=center[0],
              lon_0 = center[1]);
    # 'cyl' = equidistant cylindrical projection
    # 'aeqd' = Azimuthal Equidistant Projection              
              
fig = plt.figure()
ax = fig.add_subplot(111)

# set axis limits to basemap coordinate system
ax.set_xlim(center[1]-map_radius, center[1]+map_radius)
ax.set_ylim(center[0]-map_radius, center[0]+map_radius)
ax.set_aspect(1)        # square up axes

lon = [-122.2399758]
lat = [37.8058428]

patches = [PolygonPatch(Point(map(lon, lat)).buffer(1.0 * randint(7500, 15000)),
                        fc='#cc00cc', ec='#555555', alpha=0.5, zorder=4)
           for lon, lat in zip(lon, lat)]
              
map.drawmapboundary(fill_color='aqua')
map.fillcontinents(color='#ddaa66',lake_color='aqua')
map.drawcoastlines()

#map.readshapefile('../sample_files/comarques', 'comarques')

plt.show()    

#%%

