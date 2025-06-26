#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 14:53:05 2024
Based on Rawat et al., 2024
@author: ktravis1
"""
# 9 Panel plot per Prajjwal's paper https://amt.copernicus.org/preprints/amt-2024-114/
# This plot compares Pandora data before and after applying uncertainty limits and other criteria for WRMS and mx hz distance
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# Function to read Pandora data from PGN files
# 'Kind' options are 'rnvh' for SS NO2, 'rnvs' for DS NO2, 'rfuh' for SS HCHO, and 'rfus' for DS NO2
def get_Pandora(f1, kind):
  import pandas as pd
  import time
  from datetime import datetime
  from datetime import timedelta
  date_time = datetime.strptime("01/01/00 00:00:00", "%d/%m/%y %H:%M:%S")

  if kind == 'rnvh': # SS NO2
      col_names=['DateTime','FracDaysSince1100','Duration','SZA','SAZ','LZA','LAZ','PZA','PAZ','RMS','NRMS','ERMS','ERMS2',\
                 'P','T','O2Height','O2O2Height','O2Surf','O2O2Surf','O2Col','O2O2Col','Type','Vers','Valid','MeanVal','WTEff',\
                     'StrayLight','ShiftL1','ShiftTotal','ResChange','IntTime','BrightCount','FW1','FW2','AtmosVar','L1QF','DQ1F',\
                         'DQ2F','L2FitFlag','L2FitDQ1QF','L2FitDQ2QF','L2QFH2O','L2DQ1QFH2O','L2DQ2QFH2O','H2Osurf','IUncertH2O','H2OSurfI',\
                             'H2OHet','H2Ocol','H2OColUncert','HdistH2O','VdistH2O','L2QF','L2Q1QF','L2Q2QF','NO2surf','IUncertSurf','Mixed',\
                                 'Het','NO2strat','NO2StratUncert','NO2Col','IUncert','Hdist','Vdist',
                                 'TopHeightH2OL1','PartialH2OL1','TopHeight1','Layer1'	,'TopHeightH2OL2', 'PartialH2OL2', 'TopHeight2','Layer2', 'TopHeightH2OL3', 'PartialH2OL3', 'TopHeight3', 'Layer3', 'TopHeightH2OL4', 'PartialH2OL4', 'TopHeight4', 'Layer4', 'TopHeightH2OL5', 'PartialH2OL5', 'TopHeight5', 'Layer5', 'TopHeightH2OL6', 'PartialH2OL6', 'TopHeight6', 'Layer6', 'TopHeightH2OL7', 'PartialH2OL7', 'TopHeight7', 'Layer7', 'TopHeightH2OL8', 'PartialH2OL8', 'TopHeight8', 'Layer8', 'TopHeightH2OL9', 'PartialH2OL9', 'TopHeight9', 'Layer9', 'TopHeightH2OL10', 'PartialH2OL10', 'TopHeight10', 'Layer10', 'TopHeightH2OL11', 'PartialH2OL11', 'TopHeight11', 'Layer11', 'TopHeightH2OL12', 'PartialH2OL12', 'TopHeight12', 'Layer12', 'TopHeightH2OL13', 'PartialH2OL13', 'TopHeight13', 'Layer13']
      pandata = pd.read_csv(f1, encoding='unicode_escape', skiprows=94, delimiter='\t', names=col_names)         
      pandata['NO2Col_DU'] =pandata['NO2Col']*6.022E23/1E4/2.69E16 # moles/m2 to DU     
      pandata['NO2Col_DU'] =pandata['NO2Col']*6.022E23/1E4/2.69E16 # moles/m2 to DU     
  if kind == 'rnvs': # DS NO2
      col_names=["DateTime","FracDaysSince1100","Duration","SZA","SAZ","LZA","LAZ","RMS","NRMS","ERMS","ERMS2","P","Type","Cal","CalValid","MeanFit","WTEff","StrayLight","ShiftL1","ShiftTotal","ResChange","IntTime","BrightCount","FW1","FW2","AtmosVar","AODs","AODc","AODe","L1QF","DQ1F","DQ2F","L2FitFlag","L2FitQ1QF","L2FitQ2QF","L2QF","L2Q1QF","L2Q2QF","NO2Col","IUncert","SUncert","Cuncert","Tuncert","RMSUncert","NO2Teff","Teff1","Teff2","Teff3","Teff4","AMF","AMFUncert","DiffuseCorr","Strat","StartUncert"]
      pandata = pd.read_csv(f1, encoding='unicode_escape', skiprows=93, delimiter=' ', names=col_names)          
      pandata['NO2Col_DU'] =pandata['NO2Col']*6.022E23/1E4/2.69E16 # moles/m2 to DU     

  if kind == 'rfuh': # SS HCHO
      col_names=['DateTime','FracDaysSince1100','Duration','SZA','SAZ','LZA','LAZ','PZA','PAZ','RMS','NRMS','ERMS','ERMS2','P','T','O2Height','O2O2Height','O2Surf','O2O2Surf','O2Col','O2O2Col','Type','Vers','Valid','MeanVal','WTEff','StrayLight','ShiftL1','ShiftTotal','ResChange','IntTime','BrightCount','FW1','FW2','AtmosVar','L1QF','DQ1F','DQ2F','L2FitFlag','L2FitDQ1QF','L2FitDQ2QF','L2QF','L2DQ1QF','L2DQ2QF','HCHOsurf','IUncert','Mixed','Het','HCHOCol','HCHOColUncert','Hdist','Vdist','TopHeight1','Layer1','TopHeight2','Layer2','TopHeight3','Layer3','TopHeight4','Layer4','TopHeight5','Layer5','TopHeight6','Layer6','TopHeight7','Layer7','TopHeight8','Layer8','TopHeight9','Layer9','TopHeight10','Layer10','TopHeight11','Layer11','TopHeight12','Layer12','TopHeight13','Layer13','Topheight14','LayerHeight14']
      pandata = pd.read_csv(f1, encoding='unicode_escape', skiprows=78, delimiter=' ', names=col_names)     
      pandata['HCHOCol_DU'] =pandata['HCHOCol']*6.022E23/1E4/2.69E16 # moles/m2 to DU     
  if kind == 'rfus': # DS HCHO
      col_names=['DateTime'	,'FracDaysSince1100','Duration','SZA','SAZ','LZA','LAZ','RMS','NRMS','ERMS','ERMS2','P','Type','Vers','Valid','MeanVal','WTEff','StrayLight','ShiftL1','ShiftTotal','ResChange','IntTime','BrightCount','FW1','FW2','AtmosVar','AOD1','AOD2','AOD3','L1QF','DQ1F','DQ2F','L2FitFlag','L2FitDQ1QF','L2FitDQ2QF','L2QF','L2DQ1QF','L2DQ2QF','HCHOCol','IUncert','StructUncert','CommonUncert','HCHOColUncert','RMSUncert','Teff','TeffUnc','TeffUnc2','TeffUnc3','TUncert','AMF','AMFuncert','DiffCorr']
      pandata = pd.read_csv(f1, encoding='unicode_escape', skiprows=78, delimiter=' ', names=col_names)     
      pandata['HCHOCol_DU'] =pandata['HCHOCol']*6.022E23/1E4/2.69E16 # moles/m2 to DU     
  
  pandata['IUncert_DU'] =pandata['IUncert']*6.022E23/1E4/2.69E16 # moles/m2 to DU     

  pandata["Date"] = date_time
  pandata["Year"] = 99999
  pandata["Month"] = 99999

  dd = pandata['FracDaysSince1100']
  for d in range(len(dd)):
      tmp = date_time+ timedelta(days=(dd[d]))
      pandata["Date"][d] =tmp
      pandata['Year'][d] = pandata['Date'][d].year 
      pandata['Month'][d] = pandata['Date'][d].month
     
  return pandata


# Analyze Manhattan HCHO
f1='/Users/ktravis1/OneDrive - NASA/Proposals/TEMPO_2024/Data/Pandora135s1_ManhattanNY-CCNY_L2_rfuh5p1-8.txt'
f2='/Users/ktravis1/OneDrive - NASA/Proposals/TEMPO_2024/Data/Pandora135s1_ManhattanNY-CCNY_L2_rfus5p1-8.txt'


# -------- Get SS and DS data
# SS
pandoraSS = get_Pandora(f1,'rfuh')
pandoraSS['PctUncert'] = pandoraSS['IUncert_DU']/pandoraSS['HCHOCol_DU']
# Direct Sun
pandoraDS = get_Pandora(f2,'rfus')
pandoraDS['PctUncert'] = pandoraDS['IUncert_DU']/pandoraDS['HCHOCol_DU']

pandoraDS['Date'].min()
pandoraDS['Date'].max()

pandoraDS['Year'] = pandoraDS['Date'].dt.year
pandoraSS['Year'] = pandoraSS['Date'].dt.year

# STEP 1: mean + 3 sigma of the high quality data 
# Let's just look at 2023 since I am interested in AEROMMA/STAQS field campaign

# STEP 1a: Get high quality data over period of interest. 
pandoraDS_HQ = pandoraDS[(pandoraDS['L2QF'] == 10 ) & ( pandoraDS['IUncert_DU'] >= 0) & (pandoraDS['Year'] == 2023)]
# Added last criteria to remove weird outpiers of IUNCERT_DU < 5 
pandoraSS_HQ = pandoraSS[(pandoraSS['L2QF'] == 10 ) & ( pandoraSS['IUncert_DU'] >= 0) & ( pandoraSS['IUncert_DU'] < 5 ) & 
                         (pandoraSS['Year'] == 2023)]

# Step 1b: Get the mean + 3sigma of the independent uncertainty
filterDS = pandoraDS_HQ['IUncert_DU'].mean() + 3*pandoraDS_HQ['IUncert_DU'].std()
filterMD = pandoraSS_HQ['IUncert_DU'].mean() + 3*pandoraSS_HQ['IUncert_DU'].std()

#Column 40: Independent uncertainty of formaldehyde total vertical column amount [moles per square meter],
# -1=cross section is zero in this wavelength range, -3=spectral fitting was done, 
#but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given,
# -9=spectral fitting was not successful
pandoraDS = pandoraDS.sort_values(by='Date')
pandoraSS = pandoraSS.sort_values(by='Date')

# Step 2: Merge the original pandora DS and SS data within 5 minutes
merged_df = pd.merge_asof(pandoraDS, pandoraSS, on='Date', tolerance=pd.Timedelta('5min'), \
                          direction='nearest', suffixes=('_MD', '_DS'))
merged_df['Hour'] = merged_df['Date'].dt.hour

# Just look at the summer which is my time period of interest
merged_df_ascent = merged_df[(merged_df['Date'] > '2023-06-01') & (merged_df['Date'] < '2023-10-01')]
# Define the combinations of quality flags for filtering
combinations = [
    (10,10), (10, 11), (10, 12),
    (11, 10), (11, 11), (11, 12),
    (12, 10), (12, 11), (12, 12)
]

# Step 3: now make the 9 panel plot with the raw data
fig, axes = plt.subplots(3, 3, figsize=(15, 15))

# Plot each combination
for ax, (ds_val, md_val) in zip(axes.flatten(), combinations):
    filtered_df = merged_df_ascent[(merged_df_ascent['L2QF_DS'] == ds_val) & (merged_df_ascent['L2QF_MD'] == md_val)]
    ax.scatter(filtered_df['HCHOCol_DU_DS'], filtered_df['HCHOCol_DU_MD'])
    # Calculate and annotate correlation coefficient
    if not filtered_df.empty:
       r_value  = np.corrcoef(filtered_df['HCHOCol_DU_DS'], filtered_df['HCHOCol_DU_MD'])[0, 1]
       num_points = len(filtered_df)

#       ax.text(0.05, 0.95, f'\n$R^2$={r_value**2:.2f}', transform=ax.transAxes, fontsize=22, verticalalignment='top')
       ax.text(0.05, 0.95, f'Points: {num_points}\n$R^2$={r_value**2:.2f}', transform=ax.transAxes, fontsize=12, verticalalignment='top')

    ax.set_title(f'L2QF_DS = {ds_val}, L2QF_MD = {md_val}', fontsize=20)
    ax.set_xlabel('HCHOCol_DS', fontsize=25)
    ax.set_ylabel('HCHOCol_MD', fontsize=25)

# Adjust layout
plt.tight_layout()

# Display the plot
plt.show()


# Step 4: now make the 9 panel plot with the filtered data
# Check to see that the correlations are good enough across all quality flags. Adjust filters if needed.
# Now use uncertainty filter, the WRMS filter, and add back in the data with <10% uncertainty.
# There didn't appear to be a need for the max horizontal distance filter. 
fig, axes = plt.subplots(3, 3, figsize=(15, 15))

# First step - apply the WRMS filter
filter1 = merged_df_ascent[(merged_df_ascent['NRMS_DS'] < 0.001) & (merged_df_ascent['NRMS_MD'] < 0.001)]
# Second step, apply the uncertainty filter
filter2 = filter1[(filter1['IUncert_DU_MD'] < filterMD) & (filter1['IUncert_DU_DS'] < filterDS)]
# Third step, find the data greater than the uncerainty filter but with < 10% total uncertainty. 
filter3 = filter1[(filter1['IUncert_DU_MD'] > filterMD) & (filter1['IUncert_DU_DS'] > filterDS) & \
                  (filter1['PctUncert_MD'] < 0.1) & (filter1['PctUncert_DS'] < 0.1)]

# We are now keeping the combination of 'filter2' and 'filter3' dataframes
combined_filter = pd.concat([filter2, filter3], axis=0, ignore_index=True)

# Plot each combination
for ax, (ds_val, md_val) in zip(axes.flatten(), combinations):
    filtered_df = combined_filter[(combined_filter['L2QF_DS'] == ds_val) & (combined_filter['L2QF_MD'] == md_val) ]
    
    
    ax.scatter(filtered_df['HCHOCol_DU_DS'], filtered_df['HCHOCol_DU_MD'])
    # Calculate and annotate correlation coefficient
    if not filtered_df.empty:
       r_value  = np.corrcoef(filtered_df['HCHOCol_DU_DS'], filtered_df['HCHOCol_DU_MD'])[0, 1]
       num_points = len(filtered_df)

#       ax.text(0.05, 0.95, f'\n$R^2$={r_value**2:.2f}', transform=ax.transAxes, fontsize=22, verticalalignment='top')
       ax.text(0.05, 0.95, f'Points: {num_points}\n$R^2$={r_value**2:.2f}', transform=ax.transAxes, fontsize=12, verticalalignment='top')

    ax.set_title(f'L2QF_DS = {ds_val}, L2QF_MD = {md_val}', fontsize=20)
    ax.set_xlabel('HCHOCol_DS', fontsize=25)
    ax.set_ylabel('HCHOCol_MD',fontsize=25)

# Adjust layout
plt.tight_layout()

# Display the plot
plt.show()
