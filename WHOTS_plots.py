# -*- coding: utf-8 -*-
"""

Trying to read an ERA5 file downloaded/subsetted with CDS API

@author: jtomf
"""


# This works with base2 env
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

'''
This could be better for the time axis:
    https://forum.marine.copernicus.eu/discussion/428/how-to-convert-netcdf-time-to-python-datetime-resolved/p1

import datetime # Python standard library datetime module
from netCDF4 import Dataset,netcdftime,num2date # http://unidata.github.io/netcdf4-python/
file_in = Dataset("file.nc","r",format="NETCDF4")
tname = "time_variable_name"
nctime = file_in.variables[tname][:] # get values
t_unit = file_in.variables[tname].units # get unit  "days since 1950-01-01T00:00:00Z"
t_cal = file_in.variables[tname].calendar
tvalue = num2date(nctime,units = t_unit,calendar = t_cal)
str_time = [i.strftime("%Y-%m-%d %H:%M") for i in tvalue] # to display dates as string
'''




plt.close("all")
__figdir__ = "WHOTS"
savefig_args = {'bbox_inches':'tight', 'pad_inches':0}

# data input
datapath="./"
# fname = 'Natl_foo.nc'
fname = 'ERA5_-158E_23N_2009.nc'

data = Dataset(fname)
data.variables


tind = 0
sst = data['sst'][tind,:,:]-273.15
atmp = data['t2m'][tind,:,:]-273.15
swh = data['swh'][tind,:,:]
U = data['u10'][tind,:,:]
V = data['v10'][tind,:,:]
time = data['time'] # hours since 1900-01-01 00:00:00.0
dtime = num2date(data['time'], data['time'].units)
# cdtime = [nc_time_axis.CalendarDateTime(item, "360_day") for item in dtime]
lon = data['longitude']
lat = data['latitude']
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
sst0 = data['sst'][:,ffy,ffx]-273.15
atmp0 = data['t2m'][:,ffy,ffx]-273.15
u0 = data['u10'][:,ffy,ffx]
v0 = data['v10'][:,ffy,ffx]

# Wave parameters are on a lower-res lat/lon grid:
ffxW = np.where(np.abs(lonW[:].data-lon_pt)==np.min(np.abs(lonW[:].data-lon_pt)))
ffyW = np.where(np.abs(latW[:].data-lat_pt)==np.min(np.abs(latW[:].data-lat_pt)))
ffxW = np.squeeze(ffxW)
ffyW = np.squeeze(ffyW)
swh0 = data['swh'][:,ffyW,ffxW]


plt.close('all')

fig = plt.figure(figsize=(6,4))
plt.subplot(5,1,1)
plt.plot(dtime[:].data,sst0)
plt.ylabel('SST ($^\circ$C)')
plt.title('Met summary at ' + str(lat[ffy]) + '$^\circ$N, ' + str(lon[ffx]) + '$^\circ$E')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,2)
plt.plot(dtime[:].data,swh0)
plt.ylabel('SWH (m)')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,3)
plt.plot(dtime[:].data,atmp0)
plt.ylabel('Air temp ($^\circ$C)')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,4)
plt.plot(dtime[:].data,np.sqrt(u0**2+v0**2))
plt.ylabel('Wind speed (m/s)')

plt.subplot(5,1,5)
plt.plot(dtime[:].data,u0)
plt.plot(dtime[:].data,v0)
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
plt.title('Air temp, '+str(dtime[tind]))
map.colorbar(mappable=None, location='right', size='5%', pad='2%', label='Air temp. ($^\circ$C)')
xpt, ypt = map(lon_pt, lat_pt)
map.plot(xpt, ypt, marker='D',color='m')
map.quiver(xpt,ypt,np.mean(u0),np.mean(v0),scale=10,scale_units='inches',color='k')

##########################################
# Load WHOTS data
uopfile=r'C:\Users\jtomf\Documents\Python\WHOTS_WG_proposal\OS_WHOTS_2009_D_M.nc'
uopdata = Dataset(uopfile)
uopdata.variables

atmpW = uopdata['AIRT'][:]
sstW = uopdata['TEMP'][:]
uW = uopdata['UWND'][:]
vW = uopdata['VWND'][:]

timeW = uopdata['TIME'] # 'days since 1950-01-01 00:00:00'
dtimeW = num2date(uopdata['TIME'], uopdata['TIME'].units)

fig = plt.figure(figsize=(6,4))
plt.subplot(5,1,1)
plt.plot(dtimeW[:].data,sstW)
plt.plot(dtime[:].data,sst0)
plt.ylabel('SST ($^\circ$C)')
plt.title('Met summary at ' + str(lat[ffy]) + '$^\circ$N, ' + str(lon[ffx]) + '$^\circ$E')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,2)
plt.plot(dtimeW[:].data,sstW)
plt.ylabel('SST ($^\circ$)')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,3)
plt.plot(dtimeW[:].data,atmpW)
plt.plot(dtime[:].data,atmp0)
plt.ylabel('Air temp ($^\circ$C)')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,4)
plt.plot(dtimeW[:].data,np.sqrt(uW**2+vW**2))
plt.plot(dtime[:].data,np.sqrt(u0**2+v0**2))
plt.ylabel('Wind speed (m/s)')

plt.subplot(5,1,5)
plt.plot(dtimeW[:].data,uW)
plt.plot(dtimeW[:].data,vW)
plt.plot(dtime[:].data,u0)
plt.plot(dtime[:].data,v0)
plt.ylabel('U, V (m/s)')


#! wget http://uop.whoi.edu/currentprojects/WHOTS/data/OS_WHOTS_2011_D_M.nc

import xarray as xr
import glob

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
#plt.plot(dtime[:].data,sst0)
plt.ylabel('SST ($^\circ$C)')
plt.title('Met summary at ' + str(lat[ffy]) + '$^\circ$N, ' + str(lon[ffx]) + '$^\circ$E')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(4,1,2)
plt.plot(timeW,dsmerged.AIRT)
#plt.plot(dtime[:].data,atmp0)
plt.ylabel('Air temp ($^\circ$C)')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(4,1,3)
plt.plot(timeW,np.sqrt(dsmerged.UWND**2+dsmerged.VWND**2))
#plt.plot(dtime[:].data,np.sqrt(u0**2+v0**2))
plt.ylabel('Wind speed (m/s)')

plt.subplot(4,1,4)
plt.plot(timeW,dsmerged.UWND)
plt.plot(timeW,dsmerged.VWND)
#plt.plot(dtime[:].data,u0)
#plt.plot(dtime[:].data,v0)
plt.ylabel('U, V (m/s)')

fig = plt.figure(figsize=(6,4))
n, bins, patches = plt.hist(np.sqrt(dsmerged.UWND**2+dsmerged.VWND**2),30,density=True,edgecolor='k') # 
plt.title('Wind speed PDF')
plt.xlabel('Wind speed')
plt.ylabel('Probablility density')
