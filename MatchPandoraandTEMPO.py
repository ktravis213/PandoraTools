#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 12:22:44 2025
Get TEMPO NO2 data near Pandora locations
Right now this code pulls from TEMPO both the tropospheric and stratospheric columns as well as the full original column.
Based on the way the retrieval is done, the sum of the tropospheric and stratospheric columns is not necessarily the original column.
@author: ktravis1
"""
import pandas as pd
import glob 
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy.spatial import cKDTree


# Assume you have a dataframe of Pandora observations
pandora_obs = pd.read_parquet('/Users/ktravis1/OneDrive - NASA/HAMAQ/nearby_pandorasDS.parquet')

# Function to read TEMPO NO2 L2 file
def read_TEMPO_NO2_L2(fn):
  '''
  Modified from read_TEMPO_HCHO_L2(fn) on the ASDC github -
  #https://oss.larc-impaqt.smce.nasa.gov/user/ktravis/lab/tree/TEMPO_formaldehyde_validation_with_Pandora_04.ipynb
  
  function read_TEMPO_NO2_L2 reads the following arrays from the
  TEMPO L2 NO2 product TEMPO_NO2_L2_V01:
    vertical_column;
    vertical_column_uncertainty;
  and returns respective fields along with coordinates of the pixels.

  If one requested variables cannot be read, all returned variables are zeroed
  '''
  var_name = 'vertical_column_troposphere'
  var_name_2 = 'vertical_column_stratosphere'
  var_unc_name = 'vertical_column_troposphere_uncertainty'
  var_QF_name = 'main_data_quality_flag'

  try:
    ds = nc.Dataset(fn)

    prod = ds.groups['product'] # this opens group product, /product, as prod

    var = prod.variables[var_name] # this reads variable vertical_column from prod (group product, /product)
    trop_NO2_column = np.array(var)
    
    var = prod.variables[var_name_2] # this reads variable vertical_column from prod (group product, /product)
    strat_NO2_column = np.array(var)
    
    fv_prod = var.getncattr('_FillValue')
    prod_unit = var.getncattr('units')

    var_unc = prod.variables[var_unc_name] # this reads variable vertical_column_uncertainty from prod (group product, /product)
    total_NO2_column_unc = np.array(var_unc)

    var_QF = prod.variables[var_QF_name] # this reads variable main_data_quality_flag from prod (group product, /product)
    total_NO2_column_QF = np.array(var_QF)
    fv_QF = var_QF.getncattr('_FillValue')

    geo = ds.groups['geolocation'] # this opens group geolocation, /geolocation, as geo

    lat = np.array(geo.variables['latitude']) # this reads variable latitude from geo (geolocation group, /geolocation) into a numpy array
    lon = np.array(geo.variables['longitude']) # this reads variable longitude from geo (geolocation group, /geolocation) into a numpy array
    fv_geo = geo.variables['latitude'].getncattr('_FillValue')
    time = np.array(geo.variables['time'] )# this reads variable longitude from geo (geolocation group, /geolocation) into a numpy array

    support = ds.groups['support_data'] # this opens group geolocation, /geolocation, as geo
    cloud = np.array(support.variables['eff_cloud_fraction'])
    
    fullcol = np.array(support.variables['vertical_column_total'])
    ds.close()

  except:
    print('variable '+var_name+' cannot be read in file '+fn)
    lat = 0.
    lon = 0.
    time = 0.
    fv_geo = 0.
    trop_NO2_column = 0.
    strat_NO2_column = 0.

    total_NO2_column_unc = 0.
    total_NO2_column_QF = 0.
    fv_prod = 0.
    fv_QF = -999
    prod_unit = ''
    cloud = 0.0
    fullcol = 0.0
    
  return lat, lon, fv_geo, time, trop_NO2_column, strat_NO2_column, fullcol,total_NO2_column_unc, total_NO2_column_QF, fv_prod, fv_QF, prod_unit, cloud

# Function to apply the cloud filter and quality flag criteria and flatten for cKDTree
def preprocess_tempo_file(tempo_lat, tempo_lon, trop_NO2,strat_NO2, full_NO2, cf, qf, time):
    # Apply cloud + QF filter once
    mask_valid = (cf < 0.2) & (qf == 0)

    # Flatten arrays AFTER applying the mask
    lat_flat = tempo_lat[mask_valid]
    lon_flat = tempo_lon[mask_valid]
    trop_no2_flat = trop_NO2[mask_valid]
    strat_no2_flat = strat_NO2[mask_valid]
    no2_flat = trop_no2_flat + strat_no2_flat # to make full column
    full_no2_flat = full_NO2[mask_valid]
    # Use the index of time along axis=0
    # Repeat the time for each pixel row to match 2D structure
    time_2d = np.repeat(time[:, np.newaxis], tempo_lon.shape[1], axis=1)
    time_flat = time_2d[mask_valid]

    coords = np.column_stack((lat_flat, lon_flat))
    # cKDTree makes a very fast set of k-dimensional points which can be used to rapidly look up the nearest neighbors of any point.
    #    https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html
    tree = cKDTree(coords)

    return tree, lat_flat, lon_flat, no2_flat, full_no2_flat, time_flat


# Function to find the TEMPO data in a 5km radius of of the Pandora location.
# I am sure there are other ways to do this!
def sample_sites_with_kdtree(tree, lat_flat, lon_flat, no2_flat, no2_full_flat, time_flat, site_locs):
    
    records = []
    site_coords = site_locs[[ 'Latitude','Longitude']].values

    for idx, (lat, lon) in enumerate(site_coords):
        station = site_locs.iloc[idx]['site_id']
        # 5 km in degrees - there might be a better way to do this
        nearby_ids = tree.query_ball_point([lat, lon], r=5 / 111)  # ~5 km in deg
        
        if nearby_ids:
            no2_vals = no2_flat[nearby_ids]
            no2_vals2 = no2_full_flat[nearby_ids]
            times = time_flat[nearby_ids]
            time_datetime = pd.to_datetime(times, unit='s', origin='1980-01-06')
            records.append({
                'station': station,
                'lat': lat,
                'lon': lon,
                'time': time_datetime.mean(),
                'TEMPO_NO2': no2_vals.mean(), # Strat + trop
                'TEMPO_NO2sd': no2_vals.std(),# Strat + trop

                'TEMPO_NO2_full':no2_vals2.mean(), # Original full column
                'TEMPO_NO2ds_full':no2_vals2.std(),# Orig full column
            })

    return records

from sklearn.neighbors import BallTree
# Could be more accurate but I haven't tried it yet
def sample_sites_with_balltree(lat_flat, lon_flat, no2_flat, no2_full_flat, time_flat, site_locs):
    coords_rad = np.deg2rad(np.column_stack((lat_flat, lon_flat)))
    tree = BallTree(coords_rad, metric='haversine')
    records = []
    for _, row in site_locs.iterrows():
        station = row['site_id']
        site_coord = np.deg2rad([[row['Latitude'], row['Longitude']]])
        inds = tree.query_radius(site_coord, r=5/6371)  # 5 km in radians
        if inds[0].size > 0:
            idx = inds[0]
            record = {
                'station': station,
                'lat': row['Latitude'],
                'lon': row['Longitude'],
                'time': pd.to_datetime(time_flat[idx], unit='s', origin='1980-01-06').mean(),
                'TEMPO_NO2': no2_flat[idx].mean(),
                'TEMPO_NO2sd': no2_flat[idx].std(),
                'TEMPO_NO2_full': no2_full_flat[idx].mean(),
                'TEMPO_NO2ds_full': no2_full_flat[idx].std()
            }
            records.append(record)
    return records


# First, get the number of unique Pandora locations in pandora_obs
site_locs = pandora_obs[['STATION(-)', 'LATITUDE(deg)', 'LONGITUDE(deg)']].drop_duplicates()
site_locs['site_id'] = site_locs['STATION(-)'].astype(str)

# List all TEMPO files
tempo_files = sorted(glob.glob('/Volumes/Apricorn_TLee/TEMPO/TEMPO_NO2_L2_V03_*'))

# Process each file
all_records = []

for file in tempo_files:
    print(f"Processing {file}")
    # First get the relevant variables from the file
    tempo_lat, tempo_lon, fv_geo, time, trop_NO2_column, strat_NO2_column, full_NO2_column, total_NO2_column_unc, total_NO2_column_QF, fv_prod, fv_QF, prod_unit, cloud= read_TEMPO_NO2_L2(file)
    # Apply the cloud filter and the quality flag and form the data into a cKDTree
    tree, lat_flat, lon_flat, no2_flat,no2_full_flat,time_flat = preprocess_tempo_file(tempo_lat, tempo_lon, trop_NO2_column, strat_NO2_column, full_NO2_column, cloud, total_NO2_column_QF, time)
    # Use the cKDTree to get the data near the Pandora site.  
    file_records = sample_sites_with_kdtree(tree, lat_flat, lon_flat, no2_flat, no2_full_flat,time_flat, site_locs)
    all_records.extend(file_records)

# Now you have the TEMPO data at each site location
df_all_sites = pd.DataFrame(all_records)



# Now merge my direct sun data within 30m

# Step 1: Make sure both time columns are datetime
df_all_sites['time'] = pd.to_datetime(df_all_sites['time'])
pandora_obs['Timestamp(UTC)'] = pd.to_datetime(pandora_obs['Timestamp(UTC)'])

# Step 2: Sort both DataFrames by time (required for merge_asof)
df_all_sites = df_all_sites.sort_values('time')
pandora_obs = pandora_obs.sort_values('Timestamp(UTC)')

# Step 3: Filter for only high quality (0 or 10) DS observations.  Could consider trying out 11 and 1.
pandora_obs_HQ = pandora_obs[pandora_obs['nitrogen_dioxide_l2_quality_flag(-)'] < 11]


# Step 4: Perform merge_asof per station, then concatenate
# Use a 30min time tolerance for merging TEMPO and Pandora
merged_list = []
for station in df_all_sites['station'].unique():
    tempo_sub = df_all_sites[df_all_sites['station'] == station]
    pandora_sub = pandora_obs_HQ[pandora_obs_HQ['STATION(-)'] == station]

    if len(pandora_sub) == 0:
        continue  # skip unmatched stations

    merged = pd.merge_asof(
        pandora_sub.sort_values('Timestamp(UTC)'),
        tempo_sub.sort_values('time'),
        left_on='Timestamp(UTC)',
        right_on='time',
        direction='nearest',
        tolerance=pd.Timedelta("30min")
    )
        
    merged_list.append(merged)

# Step 5: Combine all station merges
merged_df = pd.concat(merged_list, ignore_index=True)

# Step 6. =  Choose Pandora and TEMPO columns to average - otherwise you might have lots of DS observations for the same TEMPO scan in the 30min window.
merged_df_tempotime = merged_df.groupby('time').agg({
    'nitrogen_dioxide_vertical_column_amount(molecules/cm2)': 'mean',
    'TEMPO_NO2': 'mean',
    'TEMPO_NO2sd': 'mean',
    'station': 'first',  # or 'unique' if you want to keep the info
    'Nearest_City': 'first',
    'NOTE': 'first',
    'Distance_km': 'first',
})
