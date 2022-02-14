# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 14:29:58 2021

@author: jtomf
"""
import sys
sys.path.append('C:/Users/jtomf/Documents/Python/Tom_tools/')

# from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import datetime as dt
import Tom_tools_v1 as tt

#url = 'https://opendap.jpl.nasa.gov/opendap/OceanTemperature/modis/L3/aqua/11um/v2014.0/4km/daily/2012/275/A2012275.L3m_DAY_NSST_sst_4km.nc'
fname = 'Natl_foo.nc'

ds_sst = xr.open_dataset(fname)

##################
# from http://xarray.pydata.org/en/stable/plotting.html
# xarray plotting functionality is a thin wrapper around the popular matplotlib library. 
# Matplotlib syntax and function names were copied as much as possible, which makes for an easy 
# transition between the two. Matplotlib must be installed before xarray can plot.

SPURSlon = -(38+00.0017/60)
SPURSlat = 24+35.0247/60

fig = plt.figure(figsize=(8, 4))
ds_sst.sst.sel(latitude=slice(28,23),longitude=slice(-42,-34)).isel(time[0]).plot(cmap='coolwarm',levels=np.linspace(26,28,15))
plt.axis('tight')
plt.plot(SPURSlon, SPURSlat, 'o', color='k')
plt.title(ds_sst.time_coverage_start)

#Same as above, but zoomed in on the axes used for VIIRS NPP below
fig = plt.figure(figsize=(8, 4))
foo = ds_sst.sst.sel(lat=slice(28, 23), lon=slice(-42, -34)).plot.contourf(cmap='coolwarm', levels=np.linspace(27.15,27.85,20))
plt.axis('tight')
plt.plot(SPURSlon, SPURSlat, 'o', color='k')
plt.title(ds_sst.time_coverage_start)
plt.axis([-38.708669,  -37.26713,  24.23951,  25.3261])
plt.axis('scaled')


fig = plt.figure(figsize=(8, 4))
ds_sst.sst.sel(lat=slice(28,23),lon=slice(-42,-34)).plot.contourf(cmap='coolwarm',levels=np.linspace(26,28,15))
plt.axis('tight')
plt.plot(SPURSlon,SPURSlat,'o',color='k')
plt.title(ds_sst.time_coverage_start)

plt.savefig(__figdir__ + "/Figure1a.png", **savefig_args, dpi=600)

