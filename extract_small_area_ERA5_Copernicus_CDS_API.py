# -*- coding: utf-8 -*-
"""

Attempt to read ERA 5 data from https://cds.climate.copernicus.eu/user
https://cds.climate.copernicus.eu/api-how-to

The dataset and API code are here:
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form

Created on Wed Jan  6 18:02:24 2021

@author: jtomf
"""
'''
1/7/2021:
    Wow, this works great.  1 month of surface data at 2 hour intervals was
    about 8 GB.  I want to try subsetting in space.  To do this, go to the web
    interface and 
    '''
import cdsapi

# E-W valid range is -180, 180
lon0 = -158 # NORSE=3, WHOTS=-158
lat0 = 22.67 # NORSE=70, WHOTS=22.67
dlat = 5
dlon = 5
yr = '2011'


c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-single-levels', # DOI: 10.24381/cds.adbb2d47
    {
        'product_type': 'reanalysis',
        'variable': [
            '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
            '2m_temperature', 'mean_sea_level_pressure', 'mean_wave_direction',
            'mean_wave_period', 'sea_surface_temperature', 'significant_height_of_combined_wind_waves_and_swell',
            'surface_pressure', 'total_precipitation',
        ],
        'year': yr,
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
#        'time': [
#            '00:00', '02:00', '04:00', 
#            '06:00', '08:00', '10:00', 
#            '12:00', '14:00', '16:00', 
#            '18:00', '20:00', '22:00', 
#        ],      
        # area is N, W, S, E; valid range is 90, -180, -90, 180
        'area': [
            lat0+dlat, lon0-dlon, lat0-dlat,
            lon0+dlon,
        ],
        'format': 'netcdf',
    },
    'ERA5_'+str(round(lon0))+'E_'+str(round(lat0))+'N_'+ yr +'.nc')

'''
Contact

copernicus-support@ecmwf.int
Licence

Licence to use Copernicus Products
Publication date
2018-06-14
References

Citation

DOI: 10.24381/cds.adbb2d47
Related data

ERA5 hourly data on pressure levels from 1950 to 1978 (preliminary version)

ERA5 hourly data on pressure levels from 1979 to present

ERA5 hourly data on single levels from 1950 to 1978 (preliminary version)

ERA5 monthly averaged data on pressure levels from 1950 to 1978 (preliminary version)

ERA5 monthly averaged data on pressure levels from 1979 to present

ERA5 monthly averaged data on single levels from 1950 to 1978 (preliminary version)

ERA5 monthly averaged data on single levels from 1979 to present
'''
