import os
import sys
# import matplotlib as mpl
# use Agg backend for running without X-server
# mpl.use('Agg')
import csv
import iris
from src import data_paths
from src import diags_sea_cold_surges as csl1
from src import diags_eqwaves as eqw
from src import diags_precip_extremes as extreme


def sea_compute(varnames, control=None, expt=None, obs=None,
                cs_level1=True,
                eqw_level2=True,
                extreme_level3=True):
    '''
    This section does regional features assessment for SEAsia.
    There are 3 main sections:
    1. NDJFM mean, variance, and Cold surge types, statistics and composites
    2. Equatorial waves - Wavenumber-frequency filtering, computation of variance in each
       wave domain, ratio of wave variance to total variance, propagation composites
    3. Extreme precip statistics, % changes due to cold surge types, MJO, equatorial waves

    :param varnames: list of U850, V850, PRECIP, SST in that order
    :type varnames: list
    :param control: Baseline model
    :type control: Dictionary
    :param expt: Experiment model
    :type expt: Dictionary
    :param obs: Observations data
    :type obs: Dictionary
    :param cs_level1: statistics of cold surges, CES, ES, MS etc
    :type cs_level1: Logical
    :param eqw_level2: Equatorial waves
    :type eqw_level2: Logical
    :param extreme_level3: Extreme precip statistics
    :type extreme_level3: Logical
    :return:
    :rtype:
    '''

    print('Computation starts...')
    runs = [control, expt, obs]

    # Pick up only non-None runs
    runs = [r for r in runs if r != None]

    for run in runs:
        metrics = {}
        runid = run['runid']
        label = run['label']
        start_date = run['start_date']
        end_date = run['end_date']

        # data
        data_root = os.path.join(run['data_retrieve_dir'], runid)

        # Out Plot Directory
        out_plot_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
        if not os.path.exists(out_plot_dir):
            os.makedirs(out_plot_dir)

        if runid == 'obs':
            start_date = '2000/06/01'
            end_date = '2020/12/31'
            precip_file = os.path.join(data_root, 'GPM_PRECIP_SEA_24h_mean_2000_2020.nc')
            u850_file = os.path.join(data_root, 'ERA5_U850_SEA_24h_mean_2000_2020.nc')
            v850_file = os.path.join(data_root, 'ERA5_V850_SEA_24h_mean_2000_2020.nc')
            sst_file = os.path.join(data_root, 'NOAA_OISST_SEA_24h_mean_2000_2020.nc')

        else:
            start_date = run['start_date']
            end_date = run['end_date']
            precip_file = os.path.join(data_root, runid + '_PRECIP.pp.nc')
            u850_file = os.path.join(data_root, runid + '_U850.pp.nc')
            v850_file = os.path.join(data_root, runid + '_V850.pp.nc')
            sst_file = os.path.join(data_root, runid + '_SST.pp.nc')

        print(precip_file)
        var_cubes = []
        print(u850_file)
        if os.path.exists(u850_file):
            u_wind_850_cubes = iris.load_cube(u850_file)
            u_wind_850_cubes.long_name = 'x_wind_850hPa'
            # u_wind_850_cubes = u_wind_850_cubes.intersection(latitude=(-30, 30),
            #                                                 longitude=(0, 360))
            var_cubes.append(u_wind_850_cubes)

        if os.path.exists(v850_file):
            v_wind_850_cubes = iris.load_cube(v850_file)
            v_wind_850_cubes.long_name = 'y_wind_850hPa'
            # v_wind_850_cubes = v_wind_850_cubes.intersection(latitude=(-30, 30),
            #                                                 longitude=(0, 360))
            #print(v_wind_850_cubes)
            var_cubes.append(v_wind_850_cubes)
        print(precip_file)
        if os.path.exists(precip_file):
            precip_cubes = iris.load_cube(precip_file)
            precip_cubes.long_name = 'precipitation_flux'
            # precip_cubes = precip_cubes.intersection(latitude=(-30, 30),
            #                                         longitude=(0, 360))
            # for model data, convert them to mm/day
            if not runid == 'obs':
                precip_cubes.convert_units('kg m-2 day-1')

            # print(precip_cubes)
            var_cubes.append(precip_cubes)

        if os.path.exists(sst_file):
            sst_cubes = iris.load_cube(sst_file)
            sst_cubes.long_name = 'surface_temperature'
            # sst_cubes = sst_cubes.intersection(latitude=(-30, 30),
            #                                   longitude=(0, 360))

            # print(sst_cubes)
            #sst_cubes.convert_units('K')
            var_cubes.append(sst_cubes)

        # Cold Surge Level 1 diagnostics
        # Mean, variance, filtered variance, filt variance/total variance
        print(var_cubes)
        if cs_level1:
            cs_metrics = csl1.cold_surge_stats(u_wind_850_cubes, v_wind_850_cubes, runid=runid)
            metrics.update(cs_metrics)
            csl1.cold_surge_composites(var_cubes,
                                       cstype=['NDJF', 'CS', 'CES', 'MS', 'ES'],
                                       runid=runid)
            csl1.plot_cold_surge_composites(cstype=['NDJF', 'CS', 'CES', 'MS', 'ES'], runid=runid, label=label)

        # Level 2 diagnostics
        # Equatorial waves
        # print(precip_cubes)
        if eqw_level2:
            for cube in [precip_cubes, u_wind_850_cubes, v_wind_850_cubes]:
                print(cube)
                eqwave_metrics = eqw.eqwaves_compute(cube, out_plot_dir, runid, compute=True)
                metrics.update(eqwave_metrics)
            print(metrics)

        # Print metrics as csv
        csl1.print_dict(metrics, runid)

        if extreme_level3:
            extreme.precip_extremes_cs(precip_cubes, runid=runid, percentile=95, compute=True, plot=True)

        # write metrics to csv file for each suite_id
        with open(os.path.join(out_plot_dir, 'metrics.csv'), 'w') as fh:
            writer = csv.writer(fh)
            for metric in metrics.items():
                writer.writerow(metric)
