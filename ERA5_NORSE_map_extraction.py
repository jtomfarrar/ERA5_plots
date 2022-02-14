# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 11:02:33 2021

@author: jtomf
"""

import ERA5_extraction_tool
import numpy as np
import time
################
# This allows us to import Tom_tools_v1
import sys
sys.path.append('C:/Users/jtomf/Documents/Python/Tom_tools/')
################
import Tom_tools_v1 as tt


# N, W, S, E valid range is 90, -180, -90, 180
# could use lon0=0 lat0=0 dlon=180 dlat=90

# E-W valid range is -180, 180
lon0 = -158  # NORSE=3, WHOTS=-158
lat0 = 22.67  # NORSE=70, WHOTS=22.67
dlat = 5
dlon = 5
yrs = np.arange(2008,2010,1) # [2010]
region_name = 'WHOTS'

for yr in yrs:
    tt.tic()
    ERA5_extraction_tool.get_surface_vars(lon0, lat0, dlon, dlat, str(yr), region_name)
    tt.toc()
    time.sleep(25) # I am not sure this helps, but it seems like 
                  # rapid repeated requests may cause it to bog down
                  # For some reason, the first 2 requests are quickly processed
                  # and the third one takes very long


lon0 = -3  # NORSE=3, WHOTS=-158
lat0 = 70  # NORSE=70, WHOTS=22.67
dlat = 5
dlon = 5
region_name = 'NORSE'

for yr in yrs:
    tt.tic()
    ERA5_extraction_tool.get_surface_vars(lon0, lat0, dlon, dlat, str(yr), region_name)
    tt.toc()
    time.sleep(25) # I am not sure this helps, but it seems like 
                  # rapid repeated requests may cause it to bog down
                  # For some reason, the first 2 requests are quickly processed
                  # and the third one takes very long


lon0 = -145  # Papa
lat0 = 50  # Papa
dlat = 5
dlon = 5
region_name = 'Papa'

for yr in yrs:
    tt.tic()
    ERA5_extraction_tool.get_surface_vars(lon0, lat0, dlon, dlat, str(yr), region_name)
    tt.toc()
    time.sleep(25) # I am not sure this helps, but it seems like 
                  # rapid repeated requests may cause it to bog down
                  # For some reason, the first 2 requests are quickly processed
                  # and the third one takes very long

lon0 = -51  # NTAS 15°N, 51°W
lat0 = 15  # Papa
dlat = 5
dlon = 5
region_name = 'NTAS'

for yr in yrs:
    tt.tic()
    ERA5_extraction_tool.get_surface_vars(lon0, lat0, dlon, dlat, str(yr), region_name)
    tt.toc()
    time.sleep(25) # I am not sure this helps, but it seems like 
                  # rapid repeated requests may cause it to bog down
                  # For some reason, the first 2 requests are quickly processed
                  # and the third one takes very long

lon0 = -85  # Stratus 20°S, 85°W
lat0 = -20  # 
dlat = 5
dlon = 5
region_name = 'Stratus'

for yr in yrs:
    tt.tic()
    ERA5_extraction_tool.get_surface_vars(lon0, lat0, dlon, dlat, str(yr), region_name)
    tt.toc()
    time.sleep(25) # I am not sure this helps, but it seems like 
                  # rapid repeated requests may cause it to bog down
                  # For some reason, the first 2 requests are quickly processed
                  # and the third one takes very long
