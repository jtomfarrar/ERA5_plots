# -*- coding: utf-8 -*-
"""
Make plots of ERA5 surface conditions at Papa, WHOTS, NTAS, Stratus

Created on Fri Jan 22 19:38:04 2021

@author: jtomfarrar
jfarrar@whoi.edu
"""
import os
# I need this to make basemap work
os.environ['PROJ_LIB'] = r'C:\Users\jtomf\anaconda3\pkgs\cartopy-0.18.0-py38h2a8b5ed_8\Lib\site-packages\cartopy'
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mplt
import nc_time_axis
import xarray as xr
import glob

site_name = 'Jan_Mayan'#'Lofoten_Basin'#'Jan_Mayan'#'NORSE' #can be 'NTAS', 'WHOTS', 'Stratus', or 'Papa'

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
elif site_name=='Lofoten_Basin':
    lon_pt = 3
    lat_pt = 70
elif site_name=='Jan_Mayan':
    lon_pt = -8.1 #3
    lat_pt = 71 #70

__figdir__ = './figs/' + site_name
savefig_args = {'bbox_inches':'tight', 'pad_inches':0.2}

#Define path using the r prefix (which means raw string so that special character / should not be evaluated)
path = r"C:\Users\jtomf\Documents\Python\ERA5_plots"

# Get a list of all relevant .nc files
if site_name=='Jan_Mayan':
    filenames = glob.glob(path+'/ERA5_surface_NORSE_*.nc')
elif site_name=='Lofoten_Basin':
    filenames = glob.glob(path+'/ERA5_surface_NORSE_*.nc')
else:
    filenames = glob.glob(path+'/ERA5_surface_'+ site_name +'_*.nc')
    
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
#############################
# Met plot time serires
fig = plt.figure(figsize=(8,6))
plt.subplot(5,1,1)
plt.plot(time,sst0)
plt.ylabel('[$^\circ$C]')
plt.title('Met summary at ' +site_name + ' ('+ str(round(float(lat[ffy]),4)) + '$^\circ$N, ' + str(round(float(lon[ffx]),4)) + '$^\circ$E)')
plt.legend(['SST'])
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,2)
plt.plot(time,swh0)
plt.legend(['Signif. Wave Height'])
plt.ylabel('[m]')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,3)
plt.plot(time,atmp0)
plt.legend(['Air temp'])
plt.ylabel('[$^\circ$C]')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,4)
plt.plot(time,np.sqrt(u0**2+v0**2))
plt.legend(['Wind speed'])
plt.ylabel('[m/s]')
ax = plt.gca()
ax.set_xticklabels([])

plt.subplot(5,1,5)
plt.plot(time,u0)
plt.plot(time,v0)
plt.legend(['Zonal wind','Meridional wind'])
plt.ylabel('[m/s]')
plt.savefig(__figdir__+'_Met_summary',**savefig_args)

#############################
# Met plot time series
# Try doing plot with matplotlib API approach:
fig, axs = plt.subplots(3, 1, sharex=True)
wnd = 0
at = 1
wav = 2
axs[wav].plot(time, swh0)
axs[wav].set(ylabel='[m]')
axs[wav].legend(['Signif. Wave Height'])
axs[at].plot(time, atmp0)
axs[at].plot(time, sst0, 'r')
axs[at].plot(time, 0*atmp0, 'k--')
axs[at].set(ylabel='[$^\circ$C]')
axs[at].legend(['Air temp','SST'],loc='lower left')
axs[wnd].plot(time, np.sqrt(u0**2+v0**2))
axs[wnd].set(ylabel='[m/s]')
axs[wnd].legend(['Wind speed'])
axs[wnd].title.set_text('Met summary near ' +site_name + ' ('+ str(round(float(lat[ffy]),4)) + '$^\circ$N, ' + str(round(float(lon[ffx]),4)) + '$^\circ$E)')
plt.tight_layout()

ff=np.where((time>np.datetime64('2008-07-01')) & (time<np.datetime64('2008-12-01')))
for ax in axs:
    ax.set_xlim(time[ff[0][0]],time[ff[0][-1]])
ax.xaxis.set_major_locator(mplt.dates.MonthLocator()) 
ax.xaxis.set_minor_locator(mplt.dates.DayLocator()) 
axs[wav].set_ylim(ymax=10)
axs[at].set_ylim(ymin=-5)
axs[wnd].set_ylim(ymax=22)


fig.autofmt_xdate()
plt.savefig(__figdir__+'_Met_zoom',**savefig_args)


#############################
# Site map
lcc_params={'projection':'lcc', 'lat_1':lat_pt-5,'lat_2':lat_pt+5,'lat_0':lat_pt,'lon_0':lon_pt,'width':10*10**6,'height':10*10**6, 'resolution':'h'}
ortho_params = {'projection':'ortho','lat_0':lat_pt,'lon_0':lon_pt,'resolution':'l'}
dx=5
dy=5
cyl_params={'projection':'cyl', 'lat_1':lat_pt-5,'lat_2':lat_pt+5,'lat_0':lat_pt,'lon_0':lon_pt,'llcrnrlat':lat_pt-dy,'urcrnrlat':lat_pt+dy,'llcrnrlon':lon_pt-dx,'urcrnrlon':lon_pt+dx, 'resolution':'h'}

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
# plt.title(site_name + ' site ' + '('+ str(round(float(lat[ffy]),4)) + '$^\circ$N, ' + str(round(float(lon[ffx]),4)) + '$^\circ$E)')
#map.colorbar(mappable=None, location='right', size='5%', pad='2%', label='Air temp. ($^\circ$C)')
xpt, ypt = map(lon_pt, lat_pt)
map.plot(xpt, ypt, marker='D',color='m')
map.quiver(xpt,ypt,np.mean(u0),np.mean(v0),scale=10,scale_units='inches',color='k')
plt.savefig(__figdir__+'_map',**savefig_args)

#############################
# Wind/wave statistics
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



