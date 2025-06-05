#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 15:24:46 2025
Get Pandora data in a certain latitude/longitude box from RSIG - https://www.epa.gov/hesc/remote-sensing-information-gateway
The url originated in the RSIG tool and is modified here to adjust the time and bounding box
Written by Prajjwal Rawat, edited by ktravis1
@author: ktravis1
"""
import certifi
import ssl

from urllib.error import HTTPError
from concurrent.futures import ThreadPoolExecutor, as_completed
import urllib.request as urlrq
from datetime import datetime, timedelta
import pandas as pd
import numpy as np

# This should cover CONUS Pandora sites
lon1, lon2 = -126, -66
lat1, lat2 = 24,50
 
# Date range
start_date = datetime.strptime('2024-06-01', "%Y-%m-%d")
end_date = datetime.strptime('2024-08-31', "%Y-%m-%d")
 
# Generate list of dates
date_list = [(start_date + timedelta(days=x)).strftime("%Y-%m-%d")
              for x in range((end_date - start_date).days + 1)]
 
# The API call doesn't seem to like doing large bounding boxes + large dates. Thus the reading day by day.
# RSIG only allows you to pull one variable from Pandora at a time. If another variable is needed in addition to the column product and quality flag, it can be added in the function below.
def fetch_pandora_data(date_str, version, column_product, qf_product):
    try:
        base_url = "https://ofmpub.epa.gov/rsig/rsigserver?SERVICE=wcs&VERSION=1.0.0&REQUEST=GetCoverage"
        time_str = f"{date_str}T00:00:00Z/{date_str}T23:59:59Z"
        bbox_str = f"{lon1},{lat1},{lon2},{lat2}"
        context = ssl.create_default_context(cafile=certifi.where())

        # Column product
        url_col = f"{base_url}&COVERAGE={version}.{column_product}&TIME={time_str}&BBOX={bbox_str}&FORMAT=ascii&MINIMUM_QUALITY=low"
        resp_col = urlrq.urlopen(url_col, context=context)
        ddata = pd.read_csv(resp_col, sep='\t', encoding='unicode_escape')

        # Quality flag
        url_qf = f"{base_url}&COVERAGE={version}.{qf_product}&TIME={time_str}&BBOX={bbox_str}&FORMAT=ascii"
        resp_qf = urlrq.urlopen(url_qf, context=context)
        ddataQF = pd.read_csv(resp_qf, sep='\t', encoding='unicode_escape')

        # Merge
        merged = pd.merge(
            ddata, ddataQF,
            on=['Timestamp(UTC)', 'STATION(-)'],
            how='inner',
            suffixes=('', '_QF')
        )
        return merged

    except HTTPError:
        print(f"No data for {date_str}")
        return pd.DataFrame()
    except Exception as e:
        print(f"Error for {date_str}: {e}")
        return pd.DataFrame()

    
# ----- Direct-sun
all_dfs = []
with ThreadPoolExecutor(max_workers=8) as executor:
    futures = {
        executor.submit(fetch_pandora_data, d,
                        "pandora.L2_rnvs3p1_8",
                        "nitrogen_dioxide_vertical_column_amount",
                        "nitrogen_dioxide_l2_quality_flag"): d
        for d in date_list
    }
    for future in as_completed(futures):
        result = future.result()
        if not result.empty:
            all_dfs.append(result)

directsunNO2 = pd.concat(all_dfs, ignore_index=True)

# ----- Sky-scan
all_dfs = []
with ThreadPoolExecutor(max_workers=8) as executor:
    futures = {
        executor.submit(fetch_pandora_data, d,
                        "pandora.L2_rnvh3p1_8",
                        "tropospheric_nitrogen_dioxide",
                        "nitrogen_dioxide_l2_quality_flag"): d
        for d in date_list
    }
    for future in as_completed(futures):
        result = future.result()
        if not result.empty:
            all_dfs.append(result)

skyscanNO2 = pd.concat(all_dfs, ignore_index=True)

# Which sites did we get for direct sun?
print(directsunNO2['NOTE'].unique())


