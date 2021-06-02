#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 11:58:28 2021
Script to generate a heatmap for tern locations
@author: Noam Vogt-Vincent
"""

import numpy as np
from scipy.ndimage import gaussian_filter
import os
from datetime import datetime
import cmocean.cm as cmo
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from netCDF4 import Dataset

this_dir = os.path.dirname(os.path.realpath(__file__))

##############################################################################
# SETTINGS                                                                   #
##############################################################################
filter_rad = 1e0  # degrees
output_file = this_dir + '/bird_heatmap.nc'

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
# IMPORT BIRD DATA                                                           #
##############################################################################

# Import from csv
data_dir = this_dir + '/migration_data.csv'

# Assign unique integer IDs
id_list = np.genfromtxt(data_dir,
                        delimiter=',',
                        usecols=(1),
                        dtype='str')[1:]

old_ids = list(set(id_list))
new_ids = np.arange(0, len(old_ids), 1)
id_conv = dict(zip(old_ids, new_ids))

bird_id = np.zeros(np.shape(id_list), dtype='i2')

for i in range(len(id_list)):
    bird_id[i] = id_conv[id_list[i]]

# Convert dates to datetime format
date_list = np.genfromtxt(data_dir,
                          delimiter=',',
                          usecols=(3),
                          dtype='str')[1:]

time_list = np.genfromtxt(data_dir,
                          delimiter=',',
                          usecols=(4),
                          dtype='str')[1:]

bird_time = np.zeros(np.shape(id_list), dtype='O')
bird_time_str = np.zeros(np.shape(id_list), dtype='S')

for i in range(len(id_list)):
    time = date_list[i] + '-' + time_list[i]

    if len(time) < 19:
        time = time + ':00'

    bird_time[i] = datetime.strptime(time, '%d/%m/%Y-%H:%M:%S')
    bird_time_str[i] = time

# Calculate the time spent that point
duration = np.zeros(np.shape(id_list), dtype='f4')

for i in range(len(id_list)):
    current_id = bird_id[i]
    current_time = bird_time[i]

    if i <= 6664:
        next_id = bird_id[i+1]
        next_time = bird_time[i+1]

        if next_id == current_id:
            dt = next_time - current_time
            dt = dt.total_seconds()/86400.
        else:
            dt = 0.5  # Take the time to be 12 hours if this is last datapoint
    else:
        dt = 0.5

    duration[i] = dt


# Write to dictionary
bird = {
        'id': bird_id,
        'time': bird_time,
        'lat': np.genfromtxt(data_dir,
                             delimiter=',',
                             usecols=(5),
                             dtype='f8')[1:],
        'lon': np.genfromtxt(data_dir,
                             delimiter=',',
                             usecols=(6),
                             dtype='f8')[1:],
        'duration': duration
        }

##############################################################################
# BIN LOCATIONS                                                              #
##############################################################################

# For now, we'll just ignore time
loc_hist = np.histogram2d(bird['lon'],
                          bird['lat'],
                          [bin_lon[0, :], bin_lat[:, 0]],
                          weights = bird['duration'])[0]

# Apply the gaussian filter
loc_pdf = gaussian_filter(loc_hist, filter_rad/1., mode='wrap')


##############################################################################
# PLOTTING                                                                   #
##############################################################################

# Plot the histogram
f, a0 = plt.subplots(1, 1, figsize=(20, 10),
                     subplot_kw={'projection': ccrs.PlateCarree()})
f.subplots_adjust(hspace=0, wspace=0, top=0.925, left=0.1)

# Set up the colorbar
axpos = a0.get_position()
pos_x = axpos.x0+axpos.width + 0.22
pos_y = axpos.y0
cax_width = 0.02
cax_height = axpos.height

pos_cax = f.add_axes([pos_x, pos_y, cax_width, cax_height])

bath = a0.pcolormesh(bin_lon,
                     bin_lat,
                     loc_pdf.T,
                     cmap=cmo.tempo,
                     vmin=0,
                     vmax=2)

gl = a0.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False

cb = plt.colorbar(bath, cax=pos_cax)
cb.set_label('Tern-days', size=12)
a0.set_aspect('auto', adjustable=None)
a0.margins(x=-0.01, y=-0.01)
a0.coastlines()
a0.set_title('Heatmap of Arctic tern locations (with a 1 degree gaussian filter)')
plt.savefig(this_dir + 'heatmap.png', dpi=300)

##############################################################################
# EXPORT                                                                     #
##############################################################################

with Dataset(output_file, mode='w') as nc:
    # Create the dimensions
    nc.createDimension('lon', np.shape(lon)[1])
    nc.createDimension('lat', np.shape(lon)[0])
    nc.createDimension('data', len(bird_id))

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

    nc.createVariable('density', 'f4', ('lat', 'lon'), zlib=True)
    nc.variables['density'].long_name = 'tern_density'
    nc.variables['density'].units = 'tern-days'
    nc.variables['density'].standard_name = 'tern_density'
    nc.variables['density'].gaussian_stdev = '1 degree'
    nc.variables['density'][:] = loc_pdf.T

    nc.createVariable('bird_id', 'f4', ('data'), zlib=True)
    nc.variables['bird_id'].long_name = 'bird_id_standardised'
    nc.variables['bird_id'].units = 'no units'
    nc.variables['bird_id'].standard_name = 'bird_id'
    nc.variables['bird_id'][:] = bird['id']

    nc.createVariable('bird_longitude', 'f4', ('data'), zlib=True)
    nc.variables['bird_longitude'].long_name = 'longitude'
    nc.variables['bird_longitude'].units = 'degrees_east'
    nc.variables['bird_longitude'].standard_name = 'longitude'
    nc.variables['bird_longitude'][:] = bird['lon']

    nc.createVariable('bird_latitude', 'f4', ('data'), zlib=True)
    nc.variables['bird_latitude'].long_name = 'latitude'
    nc.variables['bird_latitude'].units = 'degrees_north'
    nc.variables['bird_latitude'].standard_name = 'latitude'
    nc.variables['bird_latitude'][:] = bird['lat']

    nc.createVariable('bird_time', 'f4', ('data'), zlib=True)
    nc.variables['bird_time'].long_name = 'time'
    nc.variables['bird_time'].units = 'no units'
    nc.variables['bird_time'].standard_name = 'time'
    nc.variables['bird_time'][:] = bird_time_str

    nc.createVariable('bird_duration', 'f4', ('data'), zlib=True)
    nc.variables['bird_duration'].long_name = 'time_spent_at_point'
    nc.variables['bird_duration'].units = 'days'
    nc.variables['bird_duration'].standard_name = 'duration'
    nc.variables['bird_duration'][:] = bird['duration']


