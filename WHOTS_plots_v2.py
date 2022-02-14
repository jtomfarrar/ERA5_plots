# -*- coding: utf-8 -*-
"""

Trying to read an ERA5 file downloaded/subsetted with CDS API

@author: jtomf
"""

import os
os.environ['PROJ_LIB'] = r'C:\Users\jtomf\anaconda3\pkgs\cartopy-0.18.0-py38h2a8b5ed_8\Lib\site-packages\cartopy'
from mpl_toolkits.basemap import Basemap

import numpy as np
from netCDF4 import Dataset
from netCDF4 import num2date
#import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib as mplt
from scipy import signal
import datetime
import nc_time_axis
import xarray as xr
import glob

plt.close("all")
__figdir__ = "WHOTS"
savefig_args = {'bbox_inches':'tight', 'pad_inches':0}

# data input
datapath="./"
# fname = 'Natl_foo.nc'

#Define path using the r prefix (which means raw string so that special character / should not be evaluated)
path_dir: str = r"C:\Users\jtomf\Documents\Python\WHOTS_WG_proposal"

# Get a list of all .nc files available in different folders
filenames = glob.glob("ERA5_surface_WHOTS_*.nc")
# filenames = [path_dir+r"\OS_WHOTS_2011_D_M.nc", path_dir+r"\OS_WHOTS_2010_D_M.nc"]
ERA = xr.open_mfdataset(filenames,combine='nested',concat_dim='time')

time = ERA.time  # 'days since 1950-01-01 00:00:00'

tind = 0
sst = ERA.sst[tind,:,:]-273.15
atmp = ERA.t2m[tind,:,:]-273.15
swh = ERA.swh[tind,:,:]
U = ERA.u10[tind,:,:]
V = ERA.v10[tind,:,:]
lon = ERA.longitude
lat = ERA.latitude
lonmesh, latmesh = np.meshgrid(lon, lat)
# There is something off with SWH-- it seems the resolution is half as fine and
#   and the data have been padded to have the same shape as SST, etc
ny, nx= np.shape(swh)
swh = swh[0:round(ny/2)+1,0:round(nx/2)+1]
lonW = lon[0:nx:2]
latW = lat[0:ny:2]
lonmeshW, latmeshW = np.meshgrid(lonW, latW)


# Pull out time series at a point
lon_pt = -158 # WHOTS=-158
lat_pt = 22.67 # WHOTS=22.67
ffx = np.where(np.abs(lon[:].data-lon_pt)==np.min(np.abs(lon[:].data-lon_pt)))
ffy = np.where(np.abs(lat[:].data-lat_pt)==np.min(np.abs(lat[:].data-lat_pt)))
ffx = np.squeeze(ffx)
ffy = np.squeeze(ffy)
sst0 = ERA.sst[:,ffy,ffx]-273.15
atmp0 = ERA.t2m[:,ffy,ffx]-273.15
u0 = ERA.u10[:,ffy,ffx]
v0 = ERA.v10[:,ffy,ffx]

# Wave parameters are on a lower-res lat/lon grid:
ffxW = np.where(np.abs(lonW[:].data-lon_pt)==np.min(np.abs(lonW[:].data-lon_pt)))
ffyW = np.where(np.abs(latW[:].data-lat_pt)==np.min(np.abs(latW[:].data-lat_pt)))
ffxW = np.squeeze(ffxW)
ffyW = np.squeeze(ffyW)
swh0 = ERA.swh[:,ffyW,ffxW]


plt.close('all')

fig = plt.figure(figsize=(6,4))
plt.subplot(5,1,1)
plt.plot(time,sst0)
plt.ylabel('SST ($^\circ$C)')
plt.title('Met summary at ' + str(round(float(lat[ffy]),4)) + '$^\circ$N, ' + str(round(float(lon[ffx]),4)) + '$^\circ$E')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,2)
plt.plot(time,swh0)
plt.ylabel('SWH (m)')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,3)
plt.plot(time,atmp0)
plt.ylabel('Air temp ($^\circ$C)')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,4)
plt.plot(time,np.sqrt(u0**2+v0**2))
plt.ylabel('Wind speed (m/s)')

plt.subplot(5,1,5)
plt.plot(time,u0)
plt.plot(time,v0)
plt.ylabel('U, V (m/s)')
# ax.set_xticklabels(dtime[:].data.strftime('%d-%m-%Y'))
# plt.xticks(rotation=-90)
# plt.title('Wind components at ' + str(lat[ffy]) + 'N, ' + str(lon[ffx]) + 'E')
plt.savefig(__figdir__+'_Met_summary',**savefig_args)

lcc_params={'projection':'lcc', 'lat_1':20.,'lat_2':30,'lat_0':lat_pt,'lon_0':lon_pt,'width':1*10**6,'height':1*10**6, 'resolution':'h'}
ortho_params = {'projection':'ortho','lat_0':35,'lon_0':lon_pt,'resolution':'l'}
dx=5
dy=5
cyl_params={'projection':'cyl', 'lat_1':20.,'lat_2':30,'lat_0':lat_pt,'lon_0':lon_pt,'llcrnrlat':lat_pt-dy,'urcrnrlat':lat_pt+dy,'llcrnrlon':lon_pt-dx,'urcrnrlon':lon_pt+dx, 'resolution':'h'}


fig = plt.figure(figsize=(12,5))
map = Basemap(**cyl_params)
x,y = map(lonmesh,latmesh) # translate lat/lon to map coordinates
map.fillcontinents()
map.drawcoastlines()
map.drawcountries()
map.contourf(x, y , atmp, cmap='coolwarm', levels=30)
# map.contourf(lonmesh, latmesh, sst, cmap='coolwarm', levels=np.linspace(-1,25,26))
map.drawparallels(range(-90, 90, 1),labels=[1,0,0,1])
map.drawmeridians(range(0, 360, 2), labels=[1,0,0,1]) #  labels=[1,0,0,1]
plt.title('Air temp, '+str(time[tind].values))
map.colorbar(mappable=None, location='right', size='5%', pad='2%', label='Air temp. ($^\circ$C)')
xpt, ypt = map(lon_pt, lat_pt)
map.plot(xpt, ypt, marker='D',color='m')
map.quiver(xpt,ypt,np.mean(u0),np.mean(v0),scale=10,scale_units='inches',color='k')

##########################################
# Load WHOTS data

#Define path using the r prefix (which means raw string so that special character / should not be evaluated)
path_dir: str = r"C:\Users\jtomf\Documents\Python\WHOTS_WG_proposal"

# Get a list of all .nc files available in different folders
filenames = glob.glob(path_dir+r"\*D_M.nc")
# filenames = [path_dir+r"\OS_WHOTS_2011_D_M.nc", path_dir+r"\OS_WHOTS_2010_D_M.nc"]
dsmerged = xr.open_mfdataset(filenames,combine='nested',concat_dim='TIME')

timeW = dsmerged.TIME  # 'days since 1950-01-01 00:00:00'

fig = plt.figure(figsize=(6,4))

plt.subplot(4,1,1)
plt.plot(timeW,dsmerged.TEMP)
plt.plot(time,sst0)
plt.ylabel('SST ($^\circ$C)')
plt.title('Met summary at ' + str(round(float(lat[ffy]),4)) + '$^\circ$N, ' + str(round(float(lon[ffx]),4)) + '$^\circ$E')
var='SST'
plt.legend(['Buoy '+var,'ERA5 '+var], loc='upper right')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(4,1,2)
plt.plot(timeW,dsmerged.AIRT)
plt.plot(time,atmp0)
plt.ylabel('Air temp ($^\circ$C)')
var='air temp'
plt.legend(['Buoy '+var,'ERA5 '+var], loc='upper right')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(4,1,3)
plt.plot(timeW,np.sqrt(dsmerged.UWND**2+dsmerged.VWND**2))
plt.plot(time,np.sqrt(u0**2+v0**2))
plt.ylabel('Wind speed (m/s)')
var='wind speed'
plt.legend(['Buoy '+var,'ERA5 '+var], loc='upper right')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(4,1,4)
plt.plot(timeW,dsmerged.UWND)
plt.plot(timeW,dsmerged.VWND)
plt.plot(time,u0)
plt.plot(time,v0)
plt.legend(['Buoy U', 'Buoy V','ERA5 U', 'ERA5 V'], loc='upper right')
plt.ylabel('U, V (m/s)')

fig = plt.figure(figsize=(6,4))
n, bins, patches = plt.hist(np.sqrt(dsmerged.UWND**2+dsmerged.VWND**2),30,density=True,edgecolor='k') # 
plt.title('Wind speed PDF')
plt.xlabel('Wind speed')
plt.ylabel('Probablility density')


# Try a climatology
dsclim = dsmerged.groupby('TIME.month', squeeze = False).mean('TIME')
fig = plt.figure(figsize=(6,4))
plt.subplot(3,1,1)
plt.plot(dsclim.TEMP)
plt.plot(dsclim.AIRT)
plt.legend(['SST','Air Temp.'], loc='upper right')
plt.ylabel('SST/Air Temp ($^\circ$C)')
plt.title('WHOTS climatology')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(3,1,2)
plt.plot(np.sqrt(dsclim.UWND**2+dsclim.VWND**2))
plt.ylabel('Wind speed (m/s)')
var='wind speed'
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(3,1,3)
plt.plot(dsclim.UWND)
plt.plot(dsclim.VWND)
plt.ylabel('U, V (m/s)')
myFmt = mplt.dates.DateFormatter('%d')
ax=plt.gca()
ax.xaxis.set_major_formatter(myFmt)

