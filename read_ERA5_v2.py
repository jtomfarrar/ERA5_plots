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
__figdir__ = "NORSE_"
savefig_args = {'bbox_inches':'tight', 'pad_inches':0}

# data input
datapath="./"
fname = 'Natl_foo.nc'
#fname = 'long_time_3E_70N.nc'

# fname = 'SPURS2_ERA5_Ray.nc'

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
lon_pt = 3
lat_pt = 70
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
plt.subplot(4,1,1)
plt.plot(dtime[:].data,sst0)
plt.ylabel('SST ($^\circ$C)')
plt.title('Met summary at ' + str(lat[ffy]) + '$^\circ$N, ' + str(lon[ffx]) + '$^\circ$E')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(4,1,2)
plt.plot(dtime[:].data,swh0)
plt.ylabel('SWH (m)')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(4,1,3)
plt.plot(dtime[:].data,atmp0)
plt.ylabel('Air temp ($^\circ$C)')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(4,1,4)
plt.plot(dtime[:].data,u0)
plt.plot(dtime[:].data,v0)
plt.ylabel('U, V (m/s)')
#ax.set_xticklabels(dtime[:].data.strftime('%d-%m-%Y'))
#plt.xticks(rotation=-90)
# plt.title('Wind components at ' + str(lat[ffy]) + 'N, ' + str(lon[ffx]) + 'E')
plt.savefig(__figdir__+'Met_summary',**savefig_args)


fig = plt.figure(figsize=(12,5))
xlim = [-10, 130]
ylim = [50, 85]
#map = Basemap(llcrnrlon=xlim[0],llcrnrlat=ylim[0],urcrnrlon=xlim[1],urcrnrlat=ylim[1], resolution = 'h', epsg=5520)
map = Basemap(projection='lcc', lat_1=45.,lat_2=55,lat_0=60,lon_0=0,width=7*10**6,height=5*10**6, resolution = 'i')
x,y = map(lonmesh,latmesh) # translate lat/lon to map coordinates
map.fillcontinents()
map.drawcoastlines()
map.drawcountries()
map.contourf(x, y , sst, cmap='coolwarm', levels=np.linspace(-1,20,30))
# map.contourf(lonmesh, latmesh, sst, cmap='coolwarm', levels=np.linspace(-1,25,26))
map.drawparallels(range(-90, 90, 10),labels=[1,0,0,1])
map.drawmeridians(range(0, 360, 10), labels=[1,0,0,1]) #  labels=[1,0,0,1]
plt.title('SST, '+str(dtime[tind]))
map.colorbar(mappable=None, location='right', size='5%', pad='2%', label='SST ($^\circ$C)')
xpt, ypt = map(lon_pt, lat_pt)
map.plot(xpt, ypt, marker='D',color='m')

fig = plt.figure(figsize=(6,5))
map = Basemap(projection='ortho',lat_0=45,lon_0=0,resolution='l')
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
#map = Basemap(llcrnrlon=xlim[0],llcrnrlat=ylim[0],urcrnrlon=xlim[1],urcrnrlat=ylim[1], resolution = 'h', epsg=5520)
x,y = map(lonmesh,latmesh) # translate lat/lon to map coordinates
map.fillcontinents()
#map.drawcoastlines()
map.contourf(x, y , sst, cmap='coolwarm', levels=np.linspace(-1,20,30))
# map.contourf(lonmesh, latmesh, sst, cmap='coolwarm', levels=np.linspace(-1,25,26))
map.drawparallels(range(-90, 90, 10),labels=[1,1,0,1])
map.drawmeridians(range(0, 360, 10), labels=[1,0,0,1]) #  labels=[1,0,0,1]
plt.title('SST, '+str(dtime[tind]))
map.colorbar(mappable=None, location='right', size='5%', pad='15%', label='SST ($^\circ$C)')
xpt, ypt = map(lon_pt, lat_pt)
map.plot(xpt, ypt, marker='D',color='m')



fig = plt.figure(figsize=(6,5))
#map = Basemap(llcrnrlon=xlim[0],llcrnrlat=ylim[0],urcrnrlon=xlim[1],urcrnrlat=ylim[1], resolution = 'h', epsg=5520)
map = Basemap(projection='lcc', lat_1=45.,lat_2=55,lat_0=60,lon_0=0,width=7*10**6,height=5*10**6, resolution = 'i')
xW,yW = map(lonmeshW,latmeshW) # translate lat/lon to map coordinates
map.fillcontinents()
map.drawcoastlines()
map.drawcountries()
map.contourf(xW, yW , swh, levels=np.linspace(0,11,30))
# map.contourf(lonmesh, latmesh, sst, cmap='coolwarm', levels=np.linspace(-1,25,26))
map.drawparallels(range(-90, 90, 10),labels=[1,1,0,0])
map.drawmeridians(range(0, 360, 10), labels=[0,0,0,1]) #  labels=[l,r,t,b]
plt.title('SWH, '+str(dtime[tind]))
map.colorbar(mappable=None, location='right', size='5%', pad='10%', label='SWH (m)')
xpt, ypt = map(lon_pt, lat_pt)
map.plot(xpt, ypt, marker='D',color='m')

fig = plt.figure(figsize=(12,5))
map = Basemap(projection='lcc', lat_1=45.,lat_2=55,lat_0=60,lon_0=0,width=7*10**6,height=5*10**6, resolution = 'i')
x,y = map(lonmesh,latmesh) # translate lat/lon to map coordinates
map.fillcontinents()
map.drawcoastlines()
map.drawcountries()
map.contourf(x, y , atmp, cmap='coolwarm', levels=30)
# map.contourf(lonmesh, latmesh, sst, cmap='coolwarm', levels=np.linspace(-1,25,26))
map.drawparallels(range(-90, 90, 10),labels=[1,0,0,1])
map.drawmeridians(range(0, 360, 10), labels=[1,0,0,1]) #  labels=[1,0,0,1]
plt.title('Air temp, '+str(dtime[tind]))
map.colorbar(mappable=None, location='right', size='5%', pad='2%', label='($^\circ$C)')
xpt, ypt = map(lon_pt, lat_pt)
map.plot(xpt, ypt, marker='D',color='m')
plt.savefig(__figdir__+'AirT_map',**savefig_args)

