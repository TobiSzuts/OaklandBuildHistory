# -*- coding: utf-8 -*-
"""
Created on Wed May 11 13:20:36 2016

@author: tobis

Reads a large geojson file and put into local memory only those parcels that are
within a certain distance of a marker.
"""

import ijson
import geopy
import geopy.distance
#import pprint
#import pandas

distance = geopy.distance.vincenty    

center = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
#center = (37.75571233182055, -122.12383758245336)       # 2nd record in file
filename = 'Alameda County Parcel Boundaries.geojson'
r = 0.5             # in km


with open(filename, 'r') as fid:
    # store complete data fields within a certain distance of tempCenter

    objects = ijson.items(fid, 'features.item')
    data = {}

#    for (i, o) in zip(range(5000), objects) :
    for (i, o) in enumerate(objects) :
        if i%1000 == 0:
            if i%10000 == 0 :
                print('{}: {} records so far'.format(i, len(data)))
            else :
                print(i, end=' ')

        pos = o['geometry']['coordinates'][0][0][0]           # get (long, lat) of first vertex
            # would be better to get centroid of the polygon
        pos = (pos[1], pos[0])  # switch to (lat, long)
        d = '{}'.format(distance(pos, center))
        d = float(d[:-3])           # convert to a number that can be compared
        if d < 0.25 :
            data[round(d,6)] = o    # 6th decimal place in lat,long is 0.1 m
#            print('Record {} is {} km away'.format(i, d))
            
print('{}: {} records so far'.format(i, len(data)))
             
tempout.close()
