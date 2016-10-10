# -*- coding: utf-8 -*-
"""
Script to:
read an input dictionary of parcel boundaries, whose key is the distance from 
a central marker,
sort by distance,
send street address to zillow for property info,
pass valid results into an output dictionary and error results into another,
remove record from original dictionary,
and save everything back to disk.

Created on Tue May 31 15:47:16 2016
@author: tobis
"""

import requests
import xmltodict
import time
import pickle

# Zillow variables keys
with open('../private/API_keys.pkl', 'rb') as datafile:
    zid = pickle.load(datafile)             # API key
zurl = 'http://www.zillow.com/webservice/GetDeepSearchResults.htm?'

inProcessFile = 'data/OaklandParcels_inProcess.pkl'
radius = 8000            # only process parcels within this radius
numToProcess = 1000      # zillow API limits to 1000 queries per day

# load data structures
with open(inProcessFile, 'rb') as fid :
    compressed = pickle.load(fid)
# ...and unpack
data_raw = compressed['data_raw']
data_queried = compressed['data_queried']
data_errors = compressed['data_errors']
del compressed

sortedKeys = [k for k in sorted(data_raw) if k < radius]

for (i, key) in zip(range(numToProcess), sortedKeys) :
    # set up timer to keep requests under 10/s
    startT = time.time()

    try :
        # read in address details from input dictionary
        zp = {'address' : '{} {} {}'.format(data_raw[key]['properties']['ADDR_HN'],
                      data_raw[key]['properties']['ADDR_SN'],
                      data_raw[key]['properties']['ADDR_ST']),
              'citystatezip' : 'Oakland, CA ' + str(data_raw[key]['properties']['ZIP']),
              'zws-id' : zid}
        r = requests.get(zurl, params=zp)
        r_dict = xmltodict.parse(r.text)['SearchResults:searchresults']

        if r_dict['message']['code'] == '0' :       # valid response?
            r_dict = r_dict['response']['results']['result']
            
            # in case the response is a list of multiple entries, take only the first one
            if type(r_dict)==list :
                r_dict = r_dict[0]
            # prune result dictionary
            r_dict.pop('links')
            r_dict.pop('zestimate')
            r_dict.pop('localRealEstate')
            data_raw[key]['zillow'] = r_dict
            # transfer  to output dictionary
            data_queried[key] = data_raw.pop(key)
            
        else :
            print('For request {}, zillow code is {}. Here''s the record:'.format(zp['address'], r_dict['message']['code']))
            print(r_dict)
            print('-'*60)
            # transfer info to error dictionary for offline analysis
            data_errors.append({'key': key, 'value': data_raw.pop(key), 
                                  'zillow': r_dict, 'source': 'zillow'})
    except Exception as exc:
        print('Unspecified error: {}'.format(exc))
        data_errors.append({'key': key, 'value': data_raw.pop(key),
                              'source': 'exception'})
            
    # log status
    print(i, ' ', zp['address'])

    endT = time.time()
    if endT - startT < 0.2 :
        time.sleep(0.2 - (endT-startT))     # rate limit to 5 calls per second

# save dictionaries back to disk
compressed = {}
compressed['data_raw'] = data_raw
compressed['data_queried'] = data_queried
compressed['data_errors'] = data_errors
with open(inProcessFile, 'wb') as fid :
    a = pickle.Pickler(fid)
    a.dump(compressed)
