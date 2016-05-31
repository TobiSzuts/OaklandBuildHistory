# -*- coding: utf-8 -*-
"""
Created on Wed May 18 10:37:32 2016

@author: tobis
"""

#%% Load 0.25 km radius data

import pickle

datafile = open('ClevelandHeights_250m_radius.pkl', 'rb')
data = pickle.load(datafile)
datafile.close()

#%% Find centroid with shapely, and recopy dictionary using this as the new key for distance

import shapely.geometry as shp
import geopy
import geopy.distance
distance = geopy.distance.vincenty    

center = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
center = (center[1], center[0])

data1 = {}
#for (i, (k, v)) in zip(range(10), data.items()) :
for (k, v) in data.items() :
    c = shp.shape(v['geometry'])
#    v['centroid'] = c.centroid.wkt
#    d = '{}'.format(distance(v['centroid'], center))       # distance can take a wkt format!
    
    v['centroid'] = (c.centroid.x, c.centroid.y)
    d = '{}'.format(distance(v['centroid'], center))
    
    d = float(d[:-3])           # convert to a number that can be compared
#    print('Old distance: {} km, new distance: {} km'.format(k, round(d,6)))
    print('Difference = {:6.2f} m'.format(abs(k-round(d,6))*1000))
    data1[round(d,6)] = v
    
data = data1
del data1


#%% How to save key to file
import pickle 

apifile = open('../private/API_keys.pkl', 'wb')
a = pickle.Pickler(apifile)
a.dump(zid)
a.dump(gid)
apifile.close()

#%% Load API keys
datafile = open('../private/API_keys.pkl', 'rb')
zid = pickle.load(datafile)
gid = pickle.load(datafile)
datafile.close()

#%% Test strategy to get addresses
# Get zillow data for control set

import pprint
import requests

Spruce_control = [635, 646, 627, 626, 621, 620,  # between McKinley and Prospect
                  592, 583, 586, 575, 580, 571, 574, 563, 568, 555]  # between Prospect and Cleveland


#Spruce_control = Spruce_control[:3]

zp = {'address' : '657 Capell Street', 'citystatezip' : 'Oakland, CA', 
     'zws-id' : zid}
zurl = 'http://www.zillow.com/webservice/GetDeepSearchResults.htm?'

z = []
for num in Spruce_control :
    zp['address'] = '{} Spruce Street'.format(num)
#    pprint.pprint(zp)

    z.append(requests.get(zurl, params=zp))
    
#%% Parse response
import xmltodict

import geopy
import geopy.distance
distance = geopy.distance.vincenty    

z2 = {}
for (add, d) in zip(Spruce_control, z) :
    z2[str(add) + ' Spruce Street'] = xmltodict.parse(d.text)['SearchResults:searchresults']
#    z2.append(minidom.parseString(x.text))     # read into a tree
#    print(minidom.parseString(x.text))
    
data_clean = {}
for (key, value) in z2.items() :
    code = value['message']['code']
    # make sure request was valid
    if code != '0' :
        print('Code was invalid: {}'.format(code))
    else :
        lat = float(value['response']['results']['result']['address']['latitude'])
        lon = float(value['response']['results']['result']['address']['longitude'])
        pos = (lat, lon)
    #    print(key, pos, sep=', ')
        
        # now iterate over parcel data to find closest match
        d = 1
        for (k2, v2) in data.items() :
            curr_pos = v2['centroid']
            curr_pos = (curr_pos[1], curr_pos[0])       # reverse order
            d_curr = distance(pos, curr_pos).km
            if d_curr < d :
                d = d_curr
                closest_parcel = v2
                closest_parcel['key'] = k2
                
        # add zillow data to master database
        data[closest_parcel['key']]['zillow'] = value['response']['results']['result']
    
        # create new database to simplify development
        data_clean[closest_parcel['key']] = data[closest_parcel['key']]
        data_clean[closest_parcel['key']]['zillow'] = value['response']['results']['result']
        
#        print(data[closest_parcel['key']])
        print('{}: parcel - zillow coordinates = {:6.4f} m'.format(key, 
              distance(data[closest_parcel['key']]['centroid'], 
                       (data[closest_parcel['key']]['zillow']['address']['longitude'],
                        data[closest_parcel['key']]['zillow']['address']['latitude'])).m))
        
        
#%% Save smaller complete database to file - this is to confirm accuracy of address
import pickle 

apifile = open('Mapping3_address.pkl', 'wb')
a = pickle.Pickler(apifile)
#a.dump(data)
a.dump(data_clean)
apifile.close()

#%% Load dev database
datafile = open('Mapping3_address.pkl', 'rb')
#data = pickle.load(datafile)
data_clean = pickle.load(datafile)
datafile.close()

#%% Query Google Maps for reverse geocode
import requests

#p = {'address' : '650 Spruce St, Oakland, CA'.replace(' ', '+')}
for (i, (key, value)) in zip(range(20), data_clean.items()) :
    

    lat = value['centroid'][1]
    long = value['centroid'][0]
    p = {'latlng' : '{},{}'.format(lat, long), 'key': gid}
    
    url = 'https://maps.googleapis.com/maps/api/geocode/json?'
    r = requests.get(url, params=p)

    # get number + full street name 
    gaddress = (r.json()['results'][0]['address_components'][0]['long_name'] + 
                ' ' + r.json()['results'][0]['address_components'][1]['long_name'])
    # save to database
    data_clean[key]['gaddress'] = gaddress
    
    print('Address = ' + value['zillow']['address']['street'] + ': Google = ' + gaddress)
    
#%%  Process parcel database to find house number and build info
"""
Loop through database, first sending geocode query to google, then taking that street
address and sending it to zillow.

Error codes that can be expected are 506 and 507
"""    

import requests
import xmltodict
import time

gurl = 'https://maps.googleapis.com/maps/api/geocode/json?'
zurl = 'http://www.zillow.com/webservice/GetDeepSearchResults.htm?'

# array to store all errors for later analysis
error_entries = [];
for (i, (key, value)) in zip(range(200), data.items()) :
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
    
#%% Save complete database to file
import pickle 

apifile = open('ClevelandHeights_250m_radius_complete.pkl', 'wb')
a = pickle.Pickler(apifile)
a.dump(data)
apifile.close()

#%% Load complete database

import pickle

datafile = open('ClevelandHeights_250m_radius_complete.pkl', 'rb')
#data = pickle.load(datafile)
data_complete = pickle.load(datafile)
datafile.close()

#%% correct errors in dataset manually

for x in error_entries :
    for (key, values) in data.items() :
        if x['key']['centroid'] == values['centroid']:
            print(key)
            
            
    
    