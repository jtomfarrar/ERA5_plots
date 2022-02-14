# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 09:56:46 2021

Test download of SPURS2 ERA5 data for Ray:
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=form

@author: jtomf
"""
import cdsapi




var=['divergence', 'specific_humidity', 'temperature']

# 'specific_cloud_liquid_water_content', 'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity'

for varname in var:
    outputfile = 'SPURS2_ERA5_Ray_' + varname + '.nc'
    c = cdsapi.Client()
    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': [varname],
            'pressure_level': [
                '1', '2', '3', '5', '7', '10',
                '20', '30', '50', '70', '100', '125',
                '150', '175', '200', '225', '250', '300',
                '350', '400', '450', '500', '550', '600',
                '650', '700', '750', '775', '800', '825',
                '850', '875', '900', '925', '950', '975',
                '1000'],
            'year': '2016',
            'time': [
                '00:00', '01:00', '02:00', '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00', '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00', '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00', '21:00', '22:00', '23:00',
            ],
            'area': [
                15, -128, 5,
                -122,
            ],
            'month': ['08', '09',],
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
        },
        outputfile)


            
            
            
            
            
