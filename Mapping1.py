# -*- coding: utf-8 -*-
"""
Created on Wed May  4 15:29:12 2016

@author: tobis

Mapping with google apis.

"""



import requests

#%%from apiclient.discovery import build
#service = build('api_name', 'api_version')
if False :
    p = {'address' : '650 Spruce St, Oakland, CA'.replace(' ', '+')}
    url = 'https://maps.googleapis.com/maps/api/geocode/json?'

    r = requests.get(url, params=p)
#    r.status_code
#    print(r.json())


#%% sample API call for Zillow
#http://www.zillow.com/webservice/GetDeepSearchResults.htm?zws-id=<ZWSID>&address=2114+Bigelow+Ave&citystatezip=Seattle%2C+WA
if False :
# zws-id should be dfined on the command line for security reasons
    zp = {'address' : '657 Capell Street', 'citystatezip' : 'Oakland, CA', 
         'zws-id' : zid}
    zurl = 'http://www.zillow.com/webservice/GetDeepSearchResults.htm?'
    
    z = requests.get(zurl, params=zp)
    #z_parse.toxml()            # to print this
    
    
    # How to extract the year built from a zillow request
    from xml.dom import minidom
    z_parse = minidom.parseString(z.text)       # read into a tree
    
    zcode = z_parse.getElementsByTagName('code')
    if zcode[0].toxml() != '<code>0</code>' :    # make sure request was valid
        print('Invalid request')
    
    year = z_parse.getElementsByTagName('yearBuilt')
    print(year[0].childNodes[0].toxml())

#%% parse Alameda County Parcel data 
if True :
    import ijson
    import geopy
    import geopy.distance
    import pprint
    import pandas
    
    distance = geopy.distance.vincenty    
    
    center = (37.8058428, -122.2399758)        # (lat, long), Armenian Church
    filename = 'Alameda County Parcel Boundaries.geojson'

    def jdefault(o) :
        # support file for treating Decimal objects as floats, needed for json.dump
        import decimal
        if isinstance(o, decimal.Decimal) :
            return float(o)
        return o.__dict__
            

    with open(filename, 'r') as fid:

        if False :
            # look at first several fields
            for (i, o) in zip(range(2), ijson.parse(fid)) :
                print(o)
        
        if False :
            # look for parcels within a certain distance of center
            tempCenter = (37.75571233182055, -122.12383758245336)
            objects = ijson.items(fid, 'features.item.geometry.coordinates.item')
            distances = []
            for (i, o) in zip(range(200), objects) :
                pos = o[0][0]           # get (long, lat) of first vertex
                pos = (pos[1], pos[0])  # switch to (lat, long)
                d = '{}'.format(distance(pos, tempCenter))
                d = float(d[:-3])           # convert to a number that can be compared
                distances.append(d)
                if d < 0.5 :
                    print('Record {} is {} km away'.format(i, d))
                    
        
        if False :
            # save complete data fields within a certain distance of tempCenter
            tempCenter = (37.75571233182055, -122.12383758245336)       # 2nd record in file
            objects = ijson.items(fid, 'features.item')
            distances = []
            tempout = pandas.HDFStore('DeleteMe.h5')
            for (i, o) in zip(range(5), objects) :
                pos = o['geometry']['coordinates'][0][0][0]           # get (long, lat) of first vertex
                pos = (pos[1], pos[0])  # switch to (lat, long)
                d = '{}'.format(distance(pos, tempCenter))
                d = float(d[:-3])           # convert to a number that can be compared
                distances.append(d)
                if d < 0.5 :
                    print('Record {} is {} km away'.format(i, d))
#                    print('-'*60)
#                    print(o)
                    # can't use distance as key becuase a float - would need to round,
                    # and then best to compute centroid
                    s = pandas.Series(o)
                    tempout[str(round(d,6))] = s         # 6th decimal place in lat,long is 0.1 m
#                    
            tempout.close()
        
        if False :
            # count how many records there are
            objects = ijson.items(fid, 'features.item.properties.fid_parcel')
            for (i, o) in enumerate(objects):
                if i%1000 == 0:
                    print(i)
                pass
            print('Total # of records', (i+1))

        if False :
            # do any parcels have street address associated with them?
            firstParcel = (37.756209609069614, -122.12036188789048)
            
            objects = ijson.items(fid, 'features.item.properties')
#            objects = ijson.items(fid, 'features.item')
            for (i, o) in enumerate(objects) :
                if i%1000 == 0:
                    if i%10000 == 0 :
                        print(i)
                    else :
                        print(i, end=' ')

                if 'address' in o:
                    print(i)
                    pprint.pprint(o)
                    break
            
            
    fid.close()    
    
    


