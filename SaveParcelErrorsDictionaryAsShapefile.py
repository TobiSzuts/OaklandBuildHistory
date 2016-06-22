# -*- coding: utf-8 -*-
"""
Save python error dictionary created by CreateParcelDatabase.py into a shapefile.
This requires that the scheme of the original shapefile is modified in the 
same way as the output dictionary.  Optionally, can filter results by distance
from central marker.

Created on Mon Jun 13 14:28:11 2016
@author: tobis
"""

import fiona
import pickle
import os
from collections import OrderedDict
import numpy as np

inProcessFile = 'data/OaklandParcels_inProcess.pkl'         # data source
with open(inProcessFile, 'rb') as fid :
    compressed = pickle.load(fid)
# ...and unpack
data_raw = compressed['data_raw']
data_queried = compressed['data_queried']
data_errors_array = compressed['data_errors']
del compressed
data_errors = {f['key']: f for f in data_errors_array}

baseFile = 'data/Oakland_parcels/parcels'       # source of shape info in dictionary
outputFileName = 'Oakland_parcels_errors'
radius = 2000                                   # in m

# create output directory if it doesn't exist yet
if os.path.isdir('data/' + outputFileName) is False :
    os.makedirs('data/' + outputFileName)
outputFile = 'data/' + outputFileName + '/' + outputFileName + '.shp'

# Register format drivers with a context manager
with fiona.drivers():
    # get schema from original file
    with fiona.open(baseFile + '.shp') as source:
        meta = source.meta
        
    # add new fields to schema file
    meta['schema']['centroid'] = ('float:19:11', 'float:19:11')
    meta['schema']['id'] = 'float:19'
    meta['schema']['type'] = 'str:50'
    meta['schema']['yearBuilt'] = 'float:10'
    meta['schema']['properties']['YEARBUILT'] = 'int:6'
    meta['schema']['properties'] = OrderedDict(meta['schema']['properties'])
    schemaOrder = meta['schema']['properties']

    with fiona.open(outputFile, 'w', **meta) as sink:
        for (i,f) in enumerate(data_errors) :
            if f <= radius :   
                data_errors[f]['properties'] = data_errors[f]['value']['properties']
                data_errors[f]['centroid'] = data_errors[f]['value']['centroid']
                data_errors[f]['id'] = data_errors[f]['value']['id']
                data_errors[f]['geometry'] = data_errors[f]['value']['geometry']
                data_errors[f]['properties']['YEARBUILT'] = np.nan
                data_errors[f].pop('key')
                try :
                    data_errors[f].pop('zillow')
                except :
                    pass
                data_errors[f].pop('value')
#                data_errors[f]['value'].pop('centroid')
                # reorder dictionary to match schema order
                data_errors[f]['properties'] = OrderedDict(
                        (k, data_errors[f]['properties'][k]) for k in schemaOrder)                
                
                sink.write(data_errors[f])
                