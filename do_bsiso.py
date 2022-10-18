import os
import sys
# import matplotlib as mpl
# use Agg backend for running without X-server
# mpl.use('Agg')
import iris
import numpy as np
import scipy.signal as signal
import pandas as pd

import src.bsiso_utils as bsiso_utils
from src import data_paths


def bsiso_compute(control=None, expt=None, obs=None,
                  season_months=[5, 6, 7, 8, 9],
                  stage1_filter_variance=True,
                  stage2_iso_peaks=True,
                  stage3_iso_lag_composite=True,
                  stage4_compute_extremes=True,
                  stage5_plot_comp=True):
    '''
    This section does the computation of native data resolution (no regridding performed)

    Computes a simple index of NH summer ISO. 5 stages to this section
    1. filters the data to 10-90? days band
    2. compute peaks of area average time series and dates
    3. lead-lag correlations with respect to peak dates
    4. Compute extremes at different phases of ISO
    5. Generate plots of composites

    :param control: Baseline model
    :type control: Dictionary
    :param expt: Experiment model
    :type expt: Dictionary
    :param obs: Observations data
    :type obs: Dictionary

    :param season_months: list of month numbers in the season
    :type season_months: List
    :param stage1_filter_variance: Compute BSISO variance
    :type stage1_filter_variance: Logical
    :param stage2_iso_peaks: Compute BSISO peak dates
    :type stage2_iso_peaks: Logical
    :param stage3_iso_lag_composite: Compute BSISO lead-lag composites
    :type stage3_iso_lag_composite: Logical
    :param stage4_compute_extremes: Compute precip extremes
    :type stage4_compute_extremes: Logical
    :param stage5_plot_comp: Generate plots
    :type stage5_plot_comp: Logical
    :return:
    :rtype:
    '''


    '''
    :param varnames: U850, V850, PRECIP, SST in that order
    :param control: Control baseline experiment
    :param expt:
    :param obs:
    :param cs_level1: Level1 statistics of cold surges, CES, ES, MS etc
    :param level2:
    :param level3:
    :return:
    '''
    print('Computation starts...')
    runs = [control, expt, obs]

    # Pick up only non-None runs
    runs = [r for r in runs if r != None]
    for run in runs:
        metrics = {}
        runid = run['runid']
        label = run['label']
        # data
        data_root = os.path.join(run['data_retrieve_dir'], runid)

        # Out Plot Directory
        out_plot_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)

        # CSV file with dates of peak ISO phases
        bsiso_peak_dates_file = os.path.join(out_plot_dir, '%s_BSISO_peak_dates.csv' % runid)
        season_dates_file = os.path.join(out_plot_dir, '%s_MJJAS_season_dates.csv' % runid)

        if not os.path.exists(out_plot_dir):
            os.makedirs(out_plot_dir)

        # For BSISO  computations, read in U850 and V850 first
        if runid == 'obs':
            start_date = '2000/06/01'
            end_date = '2020/12/31'
            precip_file = os.path.join(data_root, 'GPM_PRECIP_SEA_24h_mean_2000_2020.nc')
            olr_file = os.path.join(data_root, 'NOAA_OLR_SEA_24h_mean_2000_2020.nc')
            u850_file = os.path.join(data_root, 'ERA5_U850_SEA_24h_mean_2000_2020.nc')
            v850_file = os.path.join(data_root, 'ERA5_V850_SEA_24h_mean_2000_2020.nc')
            sst_file = os.path.join(data_root, 'NOAA_OISST_SEA_24h_mean_2000_2020.nc')
            # Vorticity and divergence to be computed
            vort850_file = os.path.join(data_root, 'ERA5_VORT850_SEA_24h_mean_2000_2020.nc')
            div850_file = os.path.join(data_root, 'ERA5_DIV850_SEA_24h_mean_2000_2020.nc')
        else:
            start_date = run['start_date']
            end_date = run['end_date']
            precip_file = os.path.join(data_root, runid + '_PRECIP.pp.nc')
            olr_file = os.path.join(data_root, runid + '_OLR.pp.nc')
            u850_file = os.path.join(data_root, runid + '_U850.pp.nc')
            v850_file = os.path.join(data_root, runid + '_V850.pp.nc')
            sst_file = os.path.join(data_root, runid + '_SST.pp.nc')
            vort850_file = os.path.join(data_root, runid + '_VORT850.pp.nc')
            div850_file = os.path.join(data_root, runid + '_DIV850.pp.nc')

        # Computing vorticity and divergence
        if not os.path.exists(vort850_file):
            u850_cube = iris.load_cube(u850_file)
            u850_cube = u850_cube.intersection(longitude=(90, 130), latitude=(-10, 25))
            u850_cube = bsiso_utils.prepare_calendar(u850_cube)

            v850_cube = iris.load_cube(v850_file)
            v850_cube = v850_cube.intersection(longitude=(90, 130), latitude=(-10, 25))
            v850_cube = bsiso_utils.prepare_calendar(v850_cube)

            bsiso_utils.compute_vort_div(u850_cube, v850_cube, vort_file=vort850_file, div_file=div850_file)
        else:
            print('Vorticity and Divergence already calculated. Skip!')

        # filter the data and compute variance and mean
        # Step 1
        if stage1_filter_variance:
            for file_name, var_name in zip([precip_file, olr_file, u850_file, v850_file, sst_file, vort850_file, div850_file],
                                           ['precipitation_flux', 'toa_outgoing_longwave_flux', 'x_wind_850', 'y_wind_850', 'sst',
                                            'vorticity_850', 'divergence_850']):
                cube = iris.load_cube(file_name)
                if var_name == 'precipitation_flux':
                    if not runid == 'obs':
                        cube.convert_units('kg m-2 day-1')

                cube = cube.intersection(longitude=(90, 130), latitude=(-10, 25))
                cube = bsiso_utils.prepare_calendar(cube)
                print(cube)

                bsiso_utils.mean_var_season(cube, varname=var_name, runid=runid,
                                            season_months=season_months)

                # Filter the data
                data_file_filt = file_name.split('.')[0] + '_filt.nc'
                if not os.path.exists(data_file_filt):
                    cube_filt = bsiso_utils.lanczos_filter(cube, window=101, low=1. / 90., high=1. / 20.)
                    iris.save(cube_filt, data_file_filt)
                    print('Written %s' % data_file_filt)
                else:
                    print('Filtered file exists. Reading it. %s' % data_file_filt)
                    cube_filt = iris.load_cube(data_file_filt)

                bsiso_utils.mean_var_season(cube_filt, varname=var_name + '_filt', runid=runid,
                                            season_months=season_months)

        # Step 2 : Iso computation
        if stage2_iso_peaks:
            # precip_filt_filename = precip_file.split('.')[0] + '_filt.nc'
            # precip_cube_filt = iris.load_cube(precip_filt_filename)
            # precip_cube_filt_dates_df = bsiso_utils.create_dates_df(precip_cube_filt)

            # redoing index with OLR
            print('Using OLR for BSISO index.')

            olr_filt_filename = olr_file.split('.')[0] + '_filt.nc'
            olr_cube_filt = iris.load_cube(olr_filt_filename)
            olr_cube_filt_dates_df = bsiso_utils.create_dates_df(olr_cube_filt)

            # Area averaged time series
            area_ts = olr_cube_filt.intersection(longitude=(100, 105), latitude=(5, 10))
            area_ts = area_ts.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)

            # Normalise the series
            ts_std = (area_ts.data - np.ma.mean(area_ts.data)) / np.ma.std(area_ts.data)

            # invert the index because it is OLR

            ts_std *= -1

            # Find peaks above 1 stdev
            peaks, _ = signal.find_peaks(ts_std, height=0.75)

            # all summer seasons
            season_df = olr_cube_filt_dates_df.loc[(olr_cube_filt_dates_df['month'].isin(season_months))]
            # Write the seasonal dates out as csv
            season_df.to_csv(season_dates_file, header=season_df.columns, index=False)
            print('Written %s' % season_dates_file)

            # select peaks in summer alone
            peaks_df = olr_cube_filt_dates_df.loc[peaks]
            peaks_df = peaks_df.loc[(peaks_df['month'].isin(season_months))]
            peaks_df.columns = ['year', 'month', 'day']

            # Write the peaks out as csv
            peaks_df.to_csv(bsiso_peak_dates_file, header=peaks_df.columns, index=False)
            print('Written %s' % bsiso_peak_dates_file)

        # Step 3 compute ISO lead-lag composite
        if stage3_iso_lag_composite:
            peaks_dates_df = pd.read_csv(bsiso_peak_dates_file)
            # Ccmposite of unfiltered data
            for file_name, var_name in zip([precip_file, olr_file, u850_file, v850_file, sst_file, vort850_file, div850_file],
                                           ['precipitation_flux', 'toa_outgoing_longwave_flux', 'x_wind_850', 'y_wind_850', 'sst',
                                            'vorticity_850', 'divergence_850']):
                cube = iris.load_cube(file_name)
                if var_name == 'precipitation_flux':
                    if not runid == 'obs':
                        cube.convert_units('kg m-2 day-1')
                bsiso_utils.make_composite(cube, peaks_dates_df, runid=runid,
                                           varname=var_name, lag=25, write_out=True)

            # Ccmposite of filtered data
            for file_name, var_name in zip([precip_file.split('.')[0] + '_filt.nc',
                                            olr_file.split('.')[0] + '_filt.nc',
                                            u850_file.split('.')[0] + '_filt.nc',
                                            v850_file.split('.')[0] + '_filt.nc',
                                            sst_file.split('.')[0] + '_filt.nc',
                                            vort850_file.split('.')[0] + '_filt.nc',
                                            div850_file.split('.')[0] + '_filt.nc'],
                                           ['precipitation_flux_filt',
                                            'toa_outgoing_longwave_flux_filt',
                                            'x_wind_850_filt',
                                            'y_wind_850_filt',
                                            'sst_filt',
                                            'vorticity_850_filt',
                                            'divergence_850_filt']):
                cube = iris.load_cube(file_name)
                bsiso_utils.make_composite(cube, peaks_dates_df, runid=runid,
                                           varname=var_name, lag=25, write_out=True)

        if stage4_compute_extremes:
            var_name = 'precipitation_flux'
            percents = [90, 95]
            cube = iris.load_cube(precip_file)

            if not runid == 'obs':
                cube.convert_units('kg m-2 day-1')

            peaks_dates_df = pd.read_csv(bsiso_peak_dates_file)
            season_dates_df = pd.read_csv(season_dates_file)

            cube_dates_df = bsiso_utils.create_dates_df(cube)

            season_inds = bsiso_utils.return_indices_of_adf(cube_dates_df, season_dates_df)
            peak_inds = bsiso_utils.return_indices_of_adf(cube_dates_df, peaks_dates_df)

            # compute the peaks composites around 3 days of the peaks
            peak_inds_around = bsiso_utils.indices_around_peaks(peak_inds, days_around_peak=3)

            for per in percents:
                peak_extreme_cube = bsiso_utils.compute_percentiles(cube, peak_inds_around, percent=per)
                season_extreme_cube = bsiso_utils.compute_percentiles(cube, season_inds, percent=per)

                # out_plot_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
                writeout_file = os.path.join(out_plot_dir,
                                             '%s_MJJAS_%s_ISO_peak_%sth_percentile.nc' % (runid, var_name, per))
                iris.save(peak_extreme_cube, writeout_file)
                print('Written %s' % writeout_file)

                writeout_file = os.path.join(out_plot_dir,
                                             '%s_MJJAS_%s_seasonal_%sth_percentile.nc' % (runid, var_name, per))
                iris.save(season_extreme_cube, writeout_file)
                print('Written %s' % writeout_file)
        # Stage 5 plot iso composites
        if stage5_plot_comp:
            bsiso_utils.plot_composite(runid, label, lat_range=(-10, 25), lon_range=(100, 120))
