#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Retrieve observations and reanalysis data \
for MJO assessment

To access the dataset one will need to create an 
account and follow instructions to install ECMWF API
https://software.ecmwf.int/wiki/display/WEBAPI/Access+MARS

"""

from ecmwfapi import ECMWFDataServer
import os
import sys
import iris
import datetime
import cf_units as unit

def retrieve_ECMWF_data(var_name, \
                  date_start='1989-01-01', date_end='2008-12-31', data_dir='/var/tmp/'):
    """Retrieve the data at 850 and 200 hPa.
    Args:
        var_name    : Variable name. E.g., x_wind 
        date_start  : string, yyyy-mm-dd format
        date_end    : string, yyyy-mm-dd format
        data_dir    : Where data to be stored
    Returns:
        data_filename
        
    Reference and examples:
        https://software.ecmwf.int/wiki/display/WEBAPI/Access+MARS
    """
    
    server = ECMWFDataServer()
    
    # var_code    : ECMWF grib code. E.g., eastward_wind = 131.128
    cf_to_ECMWF_grib_code = {'eastward_wind': 131.128,}
    
    try:
        var_code = cf_to_ECMWF_grib_code[var_name]
    except KeyError:
        print '%s key does not exist.' %var_name
    
    # ERAINT data comes as 6 hourly
    # Download them all and then generate the daily mean
    times = ['00', '06', '12', '18']
    
    # pressure level variables
    daily_filename = os.path.join(data_dir, 'obs_' + var_name + '.nc')
    
    if not os.path.exists(daily_filename):
        for time in times:
            outname = os.path.join(data_dir, '%s_%sZ.nc' % (var_name, time))
            if not os.path.exists(outname):
                server.retrieve({
                    'dataset' : "interim",
                    'step'    : "0",
                    'levtype' : "pl",
                    'date'    : "%s/to/%s" % (date_start, date_end),
                    'time'    : "%s" % time,
                    'origin'  : "all",
                    'type'    : "an",
                    'class'   : "ei",
                    'expver'  : "1",
                    'stream'  : "oper",
                    'param'   : "%s" % var_code,
                    'levelist': "850/200",
                    'area'    : "-30/0/30/357.5",
                    'grid'    : "2.5/2.5",
                    'format'  : "netcdf",
                    'target'  : "%s" % outname
                    })
        command = '/opt/ukmo/utils/bin/ncea -O ' + tmpdir + '/%s_*Z.nc ' % (var_name) + \
                   daily_filename
        print command
        res = os.system(command)
    
    if os.path.exists(daily_filename):
        # remove the 6 hourly files
        command = '/bin/rm -f ' + os.path.join(data_dir, '%s_*Z.nc ' % (var_name))
        print command
        res = os.system(command)
            
def retrieve_NOAA_OLR_data(var_name, \
                  date_start='1989-01-01', date_end='2008-12-31', data_dir='/var/tmp/'):
    """Retrieve the Outgoing Longwave Radiation (OLR) data from 
    ftp://ftp.cdc.noaa.gov/Datasets/interp_OLR/olr.day.mean.nc and 
    subset for the time range and for tropics
    Args:
        var_name    : Variable name. E.g., x_wind 
        date_start  : string, yyyy-mm-dd format
        date_end    : string, yyyy-mm-dd format
        data_dir    : Where data to be stored
    Returns:
        data_filename
        
    Reference:
        https://www.esrl.noaa.gov/psd/data/gridded/data.interp_OLR.html
        
        Liebmann B. and C.A. Smith, 1996: Description of a Complete 
        (Interpolated) Outgoing Longwave Radiation Dataset. 
        Bulletin of the American Meteorological Society, 77, 1275-1277.
         
    Citation:
        If you use the Interpolated OLR data in a publication, please 
        cite Liebmann and Smith (Bulletin of the American Meteorological 
        Society, 77, 1275-1277, June 1996).
        
        Please note: If you acquire Interpolated OLR data products from PSD, 
        we ask that you acknowledge us in your use of the data. 
        This may be done by including text such as Interpolated OLR data provided 
        by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, from their Web site at 
        http://www.esrl.noaa.gov/psd/ in any documents or publications using 
        these data. We would also appreciate receiving a copy of the relevant 
        publications. This will help PSD to justify keeping the Interpolated 
        OLR data set freely available online in the future. Thank you!
    """
    # region (Tropics)
    lon = (0, 360)
    lat = (-30, 30)
    # Extract region
    region_constraint = iris.Constraint(
            longitude=lambda cell: lon[0] <= cell <= lon[1],
            latitude=lambda cell: lat[0] <= cell <= lat[1]
            )
    
    # retrieve the data from FTP URI
    local_file = os.path.join(data_dir, 'olr.day.mean.nc') 
    if not os.path.exists(local_file):
        remote_file = 'ftp://ftp.cdc.noaa.gov/Datasets/interp_OLR/olr.day.mean.nc'
        wget_command = 'wget -c %s -O %s' % (remote_file, local_file)
        print wget_command
        res = os.system(wget_command)
        print res
    
    if os.path.exists(local_file):
        # file after cropping the data for our use
        daily_filename = os.path.join(data_dir, 'obs_' + var_name + '.nc')
        cube = iris.load_cube(local_file, 'Daily OLR')
        cube.long_name = var_name
    
        # Extract region
        cube = cube.extract(region_constraint)
        
        # Constraining time 
        # The extremely painful iris way! 
        time_units = cube.coord('time').units
        start_year, start_month, start_day = date_start.split('-')
        end_year, end_month, end_day = date_end.split('-')
        
        start_dt = datetime.datetime(int(start_year), int(start_month), int(start_day))
        start_offset = time_units.date2num(start_dt)
        
        end_dt = datetime.datetime(int(end_year), int(end_month), int(end_day))
        end_offset = time_units.date2num(end_dt)
        
        time_constraint = iris.Constraint(time=lambda cell: start_offset <= cell.point <= end_offset)  
        cube = cube.extract(time_constraint)
        
        # Save the cubes
        iris.save(cube, daily_filename, netcdf_format="NETCDF3_CLASSIC")
        
def retrieve_GPCP_precip_data(var_name, \
                  date_start='1996-10-01', date_end='2015-11-10', data_dir='/var/tmp/'):
    """Retrieve the GPCP (Daily): Global Precipitation Climatology data from 
    http://iridl.ldeo.columbia.edu/SOURCES/.NASA/.GPCP/.V1DD/.V1p2/dods
    and subset for the time range and for tropics
    Args:
        var_name    : Variable name. E.g., x_wind 
        date_start  : string, yyyy-mm-dd format
        date_end    : string, yyyy-mm-dd format
        data_dir    : Where data to be stored
    Returns:
        data_filename
        
    Reference:
        https://rda.ucar.edu/datasets/ds728.3/
        
        Huffman, G. J., R. F. Adler, M. M. Morrissey, D. T. Bolvin, S. Curtis, 
        R. Joyce, B. McGavock, and J. Susskind, 2001: Global precipitation at 
        one-degree daily resolution from multisatellite observations. 
        J. Hydrometeor., 2, 36-50 (DOI: 10.1175/1525-7541(2001)002<0036:GPAODD>2.0.CO;2).
        
        dataset_documentation: ftp://rsd.gsfc.nasa.gov/pub/1dd-v1.2/1DD_v1.2_doc.pdf
        description: GPCP 1-Degree Daily Combination (Version 1.2)
    Citation:
        Huffman, G. J., D. T. Bolvin, and R. F. Adler. 2016. GPCP Version 1.2 
        One-Degree Daily Precipitation Data Set. Research Data Archive at the 
        National Center for Atmospheric Research, Computational and Information 
        Systems Laboratory. https://doi.org/10.5065/D6D50K46.      
        Accessed* dd mmm yyyy.
        
        *Please fill in the "Accessed" date with the day, month, and year 
        (e.g. - 5 Aug 2011) you last accessed the data from the RDA.
    """
    # region (Tropics)
    lon = (0, 360)
    lat = (-30, 30)
    # Extract region
    region_constraint = iris.Constraint(
            longitude=lambda cell: lon[0] <= cell <= lon[1],
            latitude=lambda cell: lat[0] <= cell <= lat[1]
            )
    
    remote_uri = 'http://iridl.ldeo.columbia.edu/SOURCES/.NASA/.GPCP/.V1DD/.V1p2/dods'
    cube = iris.load_cube(remote_uri)
    
    daily_filename = os.path.join(data_dir, 'obs_' + var_name + '.nc')
    
    # Extract region
    cube = cube.extract(region_constraint)
    
    # Constraining time 
    # The extremely painful iris way! 
    # Data did not have a proper calendar with no units
    # so redefining a calendar
    time_units = unit.Unit('days since 1996-10-01 00:00:00', calendar=unit.CALENDAR_STANDARD)
    cube.coord('time').units = time_units    
    cube.coord('time').points = range(len(cube.coord('time').points))
    
    # Check using
    #cube.coord('time').units.num2date(0)
    start_year, start_month, start_day = date_start.split('-')
    end_year, end_month, end_day = date_end.split('-')
    
    start_dt = datetime.datetime(int(start_year), int(start_month), int(start_day))
    start_offset = time_units.date2num(start_dt)
    
    end_dt = datetime.datetime(int(end_year), int(end_month), int(end_day))
    end_offset = time_units.date2num(end_dt)
    
    time_constraint = iris.Constraint(time=lambda cell: start_offset <= cell.point <= end_offset)  
    cube = cube.extract(time_constraint)
    
    # Save the cubes
    iris.save(cube, daily_filename, netcdf_format="NETCDF3_CLASSIC")
    
if __name__ == '__main__':    
    data_dir = '/project/MJO_GCSS/hadgem3/data/obs/daily/'
    date_start = '1989-01-01'
    date_end = '2008-12-31'
    
    # 850 and 200 Pa x_wind data
    # Retrieve from ECMWF MARS system using the command below
    variable_name = 'eastward_wind'
    retrieve_ECMWF_data(variable_name, \
                        date_start=date_start, date_end=date_end, \
                        data_dir=data_dir)

    '''
    variable_name = 'toa_outgoing_longwave_flux'
    retrieve_NOAA_OLR_data(variable_name, \
                        date_start=date_start, date_end=date_end, \
                        data_dir=data_dir)
    
    variable_name = 'precipitation_flux'
    date_start = '1996-10-01'       # GPCP data starts on 1996-10-01
    date_end = '2015-11-01'       # GPCP data ends on 2015-11-01
    
    retrieve_GPCP_precip_data(variable_name, \
                        date_start=date_start, date_end=date_end, \
                        data_dir=data_dir)
    '''
