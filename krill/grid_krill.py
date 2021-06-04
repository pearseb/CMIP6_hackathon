#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 11:11:08 2021
Grid krill abundances from CPR
@author: noam
"""

import pandas as pd
import os
import numpy as np
from datetime import datetime
from netCDF4 import Dataset

this_dir = os.path.dirname(os.path.realpath(__file__))

##############################################################################
# SETTINGS                                                                   #
##############################################################################
fh1 = this_dir + '/event.txt'
fh2 = this_dir + '/occurrence.txt'
output_fh = this_dir + '/CPR_bin.nc'

##############################################################################
# IMPORT DATA                                                                #
##############################################################################

if os.path.isfile(this_dir + '/plat.npy'):
    plat = np.load(this_dir + '/plat.npy')
    plon = np.load(this_dir + '/plon.npy')
    pt = np.load(this_dir + '/pt.npy')
    pc = np.load(this_dir + '/pc.npy')
else:
    d1 = pd.read_csv(fh1, delimiter = '\t')
    d2 = pd.read_csv(fh2, delimiter = '\t')

    # Filter d2 by species (Euphausiacea)
    d2 = d2.loc[d2['taxonID'] == 88]

    # Filter by IDs also present in both
    d2 = d2.loc[d2['id'].isin(d1['id'])]
    d1 = d1.loc[d1['id'].isin(d2['id'])]

    # Sort both
    d1 = d1.sort_values(by = ['id'])
    d2 = d2.sort_values(by = ['id'])
    d1 = d1.reset_index(drop=True)
    d2 = d2.reset_index(drop=True)

    # Extract key columns
    plat = d1['decimalLatitude'].to_numpy()
    plon = d1['decimalLongitude'].to_numpy()
    pt = d1['eventDate'].to_numpy(dtype='O')

    for i in range(len(pt)):
        pt[i] = datetime.strptime(pt[i], '%Y-%m-%dT%H:%M:%SZ').year

    pt = pt.astype('i2')
    pc = d2['individualCount'].to_numpy(dtype='i2')

    d1 = []
    d2 = []

    # Convert nans to 0 in counts
    pc[np.isnan(pc)] = 0

    np.save(this_dir + '/plat.npy', plat)
    np.save(this_dir + '/plon.npy', plon)
    np.save(this_dir + '/pt.npy', pt)
    np.save(this_dir + '/pc.npy', pc)

##############################################################################
# SET UP GRID                                                                #
##############################################################################

# Set up a 1 degree grid
lon, lat = np.meshgrid(np.linspace(-179.5, 179.5, num=360, dtype='f2'),
                       np.linspace(-89.5, 89.5, num=180, dtype='f2'))

# Set up bin boundaries for the 1 degree grid
bin_lon, bin_lat = np.meshgrid(np.linspace(-180, 180, num=361, dtype='f2'),
                               np.linspace(-90, 90, num=181, dtype='f2'))


##############################################################################
# Bin counts                                                                 #
##############################################################################

yearMin = np.min(pt)
yearMax = np.max(pt)
yearCount = yearMax - yearMin + 1

hist = np.zeros((yearCount, np.shape(lon)[0], np.shape(lon)[1]))

for i in range(yearCount):
    pc_ = pc[pt == yearMin + i]
    plat_ = plat[pt == yearMin + i]
    plon_ = plon[pt == yearMin + i]
    hist[i, :, :] = np.histogram2d(plon_,
                                   plat_,
                                   [bin_lon[0, :], bin_lat[:, 0]],
                                   weights=pc_)[0].T


##############################################################################
# Export to netcdf                                                           #
##############################################################################

with Dataset(output_fh, mode='w') as nc:
    # Create the dimensions
    nc.createDimension('lon', np.shape(lon)[1])
    nc.createDimension('lat', np.shape(lon)[0])
    nc.createDimension('time', yearCount)
    nc.createDimension('data', len(plat))

    nc.createVariable('longitude', 'f4', ('lat', 'lon'), zlib=True)
    nc.variables['longitude'].long_name = 'longitude'
    nc.variables['longitude'].units = 'degrees_east'
    nc.variables['longitude'].standard_name = 'longitude'
    nc.variables['longitude'][:] = lon

    nc.createVariable('latitude', 'f4', ('lat', 'lon'), zlib=True)
    nc.variables['latitude'].long_name = 'latitude'
    nc.variables['latitude'].units = 'degrees_north'
    nc.variables['latitude'].standard_name = 'latitude'
    nc.variables['latitude'][:] = lat

    nc.createVariable('time', 'i2', ('time'), zlib=True)
    nc.variables['time'].long_name = 'time'
    nc.variables['time'].units = 'year'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'][:] = np.arange(yearMin, yearMax+1)

    nc.createVariable('count', 'i2', ('time', 'lat', 'lon'), zlib=True)
    nc.variables['count'].long_name = 'PCR_euphausiacea_count_per_bin'
    nc.variables['count'].units = 'individuals_per_bin'
    nc.variables['count'].standard_name = 'euphausiacea_count'
    nc.variables['count'][:] = hist

    nc.createVariable('PCR_lon', 'f4', ('data'), zlib=True)
    nc.variables['PCR_lon'].long_name = 'PCR_longitude'
    nc.variables['PCR_lon'].units = 'degrees_east'
    nc.variables['PCR_lon'].standard_name = 'PCR_longitude'
    nc.variables['PCR_lon'][:] = plon

    nc.createVariable('PCR_lat', 'f4', ('data'), zlib=True)
    nc.variables['PCR_lat'].long_name = 'PCR_latitude'
    nc.variables['PCR_lat'].units = 'degrees_north'
    nc.variables['PCR_lat'].standard_name = 'PCR_latitude'
    nc.variables['PCR_lat'][:] = plat

    nc.createVariable('PCR_time', 'i2', ('data'), zlib=True)
    nc.variables['PCR_time'].long_name = 'PCR_time'
    nc.variables['PCR_time'].units = 'year'
    nc.variables['PCR_time'].standard_name = 'PCR_time'
    nc.variables['PCR_time'][:] = pt

    nc.createVariable('PCR_count', 'i2', ('data'), zlib=True)
    nc.variables['PCR_count'].long_name = 'PCR_euphausiacea_count'
    nc.variables['PCR_count'].units = 'individuals'
    nc.variables['PCR_count'].standard_name = 'PCR_euphausiacea_count'
    nc.variables['PCR_count'][:] = pc
