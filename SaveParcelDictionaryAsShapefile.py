# -*- coding: utf-8 -*-
"""
Save python dictionary created by CreateParcelDatabase.py into a shapefile.
This requires that the scheme of the original shapefile is modified in the 
same way as the output dictionary.  Optionally, can filter results by distance
from central marker.

Created on Wed Jun  8 09:51:27 2016
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
del compressed

baseFile = 'data/Oakland_parcels/parcels'       # source of shape info in dictionary
outputFileName = 'Oakland_parcels_queried'
radius = 50000                                   # in m

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
        for (i,f) in enumerate(data_queried) :
            if f <= radius :   
                if 'yearBuilt' in data_queried[f]['zillow'] :
                    if float(data_queried[f]['zillow']['yearBuilt']) < 1800 :
                        print(f, data_queried[f]['zillow']['yearBuilt'])
                    data_queried[f]['properties']['YEARBUILT'] = data_queried[f]['zillow']['yearBuilt']
                    data_queried[f].pop('zillow')
                    
                else :
                    data_queried[f]['properties']['YEARBUILT'] = np.nan
                # reorder dictionary to match schema order
                data_queried[f]['properties'] = OrderedDict(
                        (k, data_queried[f]['properties'][k]) for k in schemaOrder)                
                
                sink.write(data_queried[f])
                