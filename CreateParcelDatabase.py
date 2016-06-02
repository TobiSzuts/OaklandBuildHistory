# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:47:16 2016

@author: tobis

Script to:
read a shapefile of parcel boundaries,
filter based on being in a square boundary around a marker location,
compute centroid, 
send centroid for a reverse geocode query to google,
and send street address to zillow for property info.

Error codes that can be expected are 506 and 507

This saves a shapefile to disc that has many fields in place.

WRONG APPROACH!!! 

"""

#%% Load API keys
import pickle
datafile = open('../private/API_keys.pkl', 'rb')
zid = pickle.load(datafile)
gid = pickle.load(datafile)
datafile.close()

#%% 

import requests
import xmltodict
import time
import fiona
import geopy.distance
distance = geopy.distance.vincenty    
import shapely.geometry as shp


gurl = 'https://maps.googleapis.com/maps/api/geocode/json?'
zurl = 'http://www.zillow.com/webservice/GetDeepSearchResults.htm?'


baseFile = 'data/Oakland_parcels/parcels'

center = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
radius = 0.5                          # in km
ur = distance(kilometers=radius).destination(center, +45)
ll = distance(kilometers=radius).destination(center, -135)
ur = (ur.longitude, ur.latitude)
ll = (ll.longitude, ll.latitude)
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
        meta['schema']['centroid'] = ('float', 'float')
#    sink_schema = source.schema.copy()

        # Open an output file, using the same format driver and
        # coordinate reference system as the source. The ``meta``
        # mapping fills in the keyword parameters of fiona.open().

        with fiona.open('data/Oakland_database_temp/parcels.shp', 'w', **meta) as sink:

            error_entries = []
            # Process only the records intersecting a box.
            for (f, i) in zip(source.filter(bbox=ROI), range(6)) :
                c = shp.shape(f['geometry']).centroid
                p = (c.y, c.x)

                f['centroid'] = p            
                d = distance(p, center).m
#                print('Distance = {:6.2f} m'.format(d))
#%%
                try :
                    # set parameters for google map query
                    gp = {'latlng' : '{},{}'.format(*p), 'key': gid}
                    # query
                    rg = requests.get(gurl, params=gp)
                    # query gives a lot of information on geographic context, including coordinates for city, county, etc.    
                    
                    if rg.json()['status'] == 'OK' :
                        gaddress = (rg.json()['results'][0]['address_components'][0]['long_name'] + 
                                    ' ' + rg.json()['results'][0]['address_components'][1]['short_name'])
                        # save to database
                        f['gaddress'] = gaddress
                    else :
                        print('For record {}, google status is {}. Here''s the record:'.format(f, rg.json()['status']))
                        print(rg.json())
                        print('-'*60)
                        # and log the error for analysis later
                        error_entries.append({'key': key, 'value':value,
                                              'google': rg, 'source': 'google'})
                except Exception as exc:
                    print('Unspecified error: {}'.format(exc))
                    error_entries.append({'key': value, 'source': 'exception'})
                    
            #        raise exc
                    # these variables may not exist...so ignore any errors there
                    try:        
                        error_entries[-1]['google'] = rg.json()
                        error_entries[-1]['zillow'] = r            
                    except :
                        pass
                                              
                sink.write(f)
                counter += 1



#%%



# array to store all errors for later analysis
error_entries = [];
for (i, (key, value)) in zip(range(20), data.items()) :
    startT = time.time()
    # extract centroid of this record
    lat = value['centroid'][1]
    long = value['centroid'][0]

    try :
        # set parameters for google map query
        gp = {'latlng' : '{},{}'.format(lat, long), 'key': gid}
        # query
        rg = requests.get(gurl, params=gp)
        # query gives a lot of information on geographic context, including coordinates for city, county, etc.    
        
        if rg.json()['status'] == 'OK' :
            gaddress = (rg.json()['results'][0]['address_components'][0]['long_name'] + 
                        ' ' + rg.json()['results'][0]['address_components'][1]['short_name'])
            # save to database
            value['gaddress'] = gaddress
        else :
            print('For record {}, google status is {}. Here''s the record:'.format(key, rg.json()['status']))
            print(rg.json())
            print('-'*60)
            # and log the error for analysis later
            error_entries.append({'key': key, 'value':value,
                                  'google': rg, 'source': 'google'})
            

        zp = {'address' : gaddress, 'citystatezip' : 'Oakland, CA', 
         'zws-id' : zid}
    
    #    zp['address'] = '{} Spruce Street'.format(num)
        r = requests.get(zurl, params=zp)
    
        r_dict = xmltodict.parse(r.text)['SearchResults:searchresults']
        # in case the response is a list of multiple entries, take only the first one
        if type(r_dict)==list :
            r_dict = r_dict[0]
    #    valid response?
        if r_dict['message']['code'] == '0' :
            r_dict = r_dict['response']['results']['result']
                
            r_dict.pop('links')
            r_dict.pop('zestimate')
            r_dict.pop('localRealEstate')
            value['zillow'] = r_dict
        else :
            print('For request {}, zillow code is {}. Here''s the record:'.format(gaddress, r_dict['message']['code']))
            print(r_dict)
            print('-'*60)
            error_entries.append({'key': key, 'value':value, 
                                  'google': r.json(), 'zillow': r,
                                  'source': 'zillow'})

    except Exception as exc:
        print('Unspecified error: {}'.format(exc))
        error_entries.append({'key': value, 'source': 'exception'})
        
#        raise exc
        # these variables may not exist...so ignore any errors there
        try:        
            error_entries[-1]['google'] = rg.json()
            error_entries[-1]['zillow'] = r            
        except :
            pass
            
    print(i, ' ', gaddress)
    
    endT = time.time()
    if endT - startT < 0.2 :
        time.sleep(0.2 - (endT-startT))     # rate limit to 5 calls per second, half as fast as possible 
