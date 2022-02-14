# -*- coding: utf-8 -*-
"""
Make plots of ERA5 surface conditions at Papa, WHOTS, NTAS, Stratus

Created on Fri Jan 22 19:38:04 2021

@author: jtomfarrar
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

site_name = 'Stratus' #'NTAS'

if site_name=='WHOTS':
    lon_pt = -158 # WHOTS=-158
    lat_pt = 22.67 # WHOTS=22.67
elif site_name=='Papa':
    lon_pt = -145  # Papa
    lat_pt = 50  # Papa
elif site_name=='NTAS':
    lon_pt = -51  # NTAS 15°N, 51°W
    lat_pt = 15  # NTAS
elif site_name=='Stratus':
    lon_pt = -85  # Stratus 2_pt°S, 85°W
    lat_pt = -20  # 

plt.close("all")
__figdir__ = site_name
savefig_args = {'bbox_inches':'tight', 'pad_inches':0}

# data input
datapath="./"
# fname = 'Natl_foo.nc'

#Define path using the r prefix (which means raw string so that special character / should not be evaluated)
path_dir: str = r"C:\Users\jtomf\Documents\Python\WHOTS_WG_proposal"

# Get a list of all .nc files available in different folders
filenames = glob.glob('ERA5_surface_'+ site_name +'_*.nc')
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
plt.title('Met summary at ' +site_name + ' ('+ str(round(float(lat[ffy]),4)) + '$^\circ$N, ' + str(round(float(lon[ffx]),4)) + '$^\circ$E)')
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

lcc_params={'projection':'lcc', 'lat_1':20.,'lat_2':30,'lat_0':lat_pt,'lon_0':lon_pt,'width':10*10**6,'height':10*10**6, 'resolution':'h'}
ortho_params = {'projection':'ortho','lat_0':35,'lon_0':lon_pt,'resolution':'l'}
dx=5
dy=5
cyl_params={'projection':'cyl', 'lat_1':20.,'lat_2':30,'lat_0':lat_pt,'lon_0':lon_pt,'llcrnrlat':lat_pt-dy,'urcrnrlat':lat_pt+dy,'llcrnrlon':lon_pt-dx,'urcrnrlon':lon_pt+dx, 'resolution':'h'}


fig = plt.figure(figsize=(12,5))
map = Basemap(**lcc_params)
x,y = map(lonmesh,latmesh) # translate lat/lon to map coordinates
map.fillcontinents(lake_color='aqua')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary(fill_color='aqua')
#map.contourf(x, y , atmp, cmap='coolwarm', levels=30)
map.drawparallels(range(-90, 90, 15),labels=[1,0,0,0]) #labels = [left,right,top,bottom]
map.drawmeridians(range(0, 360, 15), labels=[0,0,0,1]) #  labels=[1,0,0,1]
plt.title(site_name + ' site ' + '('+ str(round(float(lat[ffy]),4)) + '$^\circ$N, ' + str(round(float(lon[ffx]),4)) + '$^\circ$E)')
#map.colorbar(mappable=None, location='right', size='5%', pad='2%', label='Air temp. ($^\circ$C)')
xpt, ypt = map(lon_pt, lat_pt)
map.plot(xpt, ypt, marker='D',color='m')
map.quiver(xpt,ypt,np.mean(u0),np.mean(v0),scale=10,scale_units='inches',color='k')

plt.savefig(__figdir__+'_map',**savefig_args)


fig = plt.figure(figsize=(8,4))
plt.subplot(1,2,1)
n, bins, patches = plt.hist(np.sqrt(u0**2+v0**2),30,density=True,edgecolor='k') # 
plt.title('ERA5 ' + site_name + '('+ str(round(float(lat[ffy]),4)) + '$^\circ$N, ' + str(round(float(lon[ffx]),4)) + '$^\circ$E)' )
plt.xlabel('Wind speed (m/s)')
plt.ylabel('Probablility density')

plt.subplot(1,2,2)
n, bins, patches = plt.hist(swh0,30,density=True,edgecolor='k') # 
# plt.title('ERA5 ' + site_name + '('+ str(round(float(lat[ffy]),4)) + '$^\circ$N, ' + str(round(float(lon[ffx]),4)) + '$^\circ$E)' + ' SWH PDF')
plt.xlabel('Signif. Wave Height (m)')
plt.ylabel('Probablility density')
plt.savefig(__figdir__+'_stats',**savefig_args)



