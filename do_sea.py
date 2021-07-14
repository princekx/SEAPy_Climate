import os
import sys
#import matplotlib as mpl
# use Agg backend for running without X-server
#mpl.use('Agg')

import iris
from src import data_paths
from src import diags_sea_cold_surges as csl1
from src import diags_eqwaves as eqw
from src import diags_precip_extremes as extreme

def sea_compute(varnames, control=None, expt=None, obs=None,
                cs_level1=True,
                eqw_level2=True,
                mjo_level3=False,
                extreme_level4=True):
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
        start_date = run['start_date']
        end_date = run['end_date']

        # data
        data_root = os.path.join(run['data_retrieve_dir'], runid)

        # Out Plot Directory
        out_plot_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
        if not os.path.exists(out_plot_dir):
            os.makedirs(out_plot_dir)

        # For Cold surge computations, read in U850 and V850 first
        u850_file = os.path.join(data_root, runid + '_U850.pp.nc')
        u200_file = os.path.join(data_root, runid + '_U200.pp.nc')
        v850_file = os.path.join(data_root, runid + '_V850.pp.nc')
        precip_file = os.path.join(data_root, runid + '_PRECIP.pp.nc')
        sst_file = os.path.join(data_root, runid + '_SST.pp.nc')
        print(precip_file)
        var_cubes = []
        print(u850_file)
        if os.path.exists(u850_file):
            u_wind_850_cubes = iris.load_cube(u850_file)
            u_wind_850_cubes.long_name = 'x_wind_850hPa'
            u_wind_850_cubes = u_wind_850_cubes.intersection(latitude=(-30, 30),
                                                             longitude=(0, 360))
            var_cubes.append(u_wind_850_cubes)

        print(u200_file)
        if os.path.exists(u200_file):
            u_wind_200_cubes = iris.load_cube(u200_file)
            u_wind_200_cubes.long_name = 'x_wind_200hPa'
            u_wind_200_cubes = u_wind_200_cubes.intersection(latitude=(-30, 30),
                                                             longitude=(0, 360))
            print(u_wind_200_cubes)
            var_cubes.append(u_wind_200_cubes)
        else:
            print('%s not found' %u200_file)

        if os.path.exists(v850_file):
            v_wind_850_cubes = iris.load_cube(v850_file)
            v_wind_850_cubes.long_name = 'y_wind_850hPa'
            v_wind_850_cubes = v_wind_850_cubes.intersection(latitude=(-30, 30),
                                                             longitude=(0, 360))
            # print(v_wind_850_cubes)
            var_cubes.append(v_wind_850_cubes)
        print(precip_file)
        if os.path.exists(precip_file):
            precip_cubes = iris.load_cube(precip_file)
            precip_cubes.long_name = 'precipitation_flux'
            precip_cubes = precip_cubes.intersection(latitude=(-30, 30),
                                                     longitude=(0, 360))
            # for model data, convert them to mm/day
            if runid == 'obs':
                precip_cubes.convert_units('mm')
            else:
                precip_cubes.convert_units('kg m-2 day-1')

            # print(precip_cubes)
            var_cubes.append(precip_cubes)

        if os.path.exists(sst_file):
            sst_cubes = iris.load_cube(sst_file)
            sst_cubes.long_name = 'surface_temperature'
            sst_cubes = sst_cubes.intersection(latitude=(-30, 30),
                                               longitude=(0, 360))

            # print(sst_cubes)
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
            csl1.plot_cold_surge_composites(cstype=['NDJF', 'CS', 'CES', 'MS', 'ES'], runid=runid)

        # Level 2 diagnostics
        # Equatorial waves
        print(precip_cubes)
        if eqw_level2:
            for cube in [precip_cubes, u_wind_850_cubes, v_wind_850_cubes]:
                print(cube)
                eqwave_metrics = eqw.eqwaves_compute(cube, out_plot_dir, runid, compute=True)
                metrics.update(eqwave_metrics)
            print(metrics)

        # Print metrics as csv
        csl1.print_dict(metrics, runid)
        '''
        # Level 3 diagnostics
        # Real-time multivariate MJO Index (RMM) calculations, and Summer/Winter
        # composites
        if level3:
            level3_metrics = diags_level3.diagnos_level3(
                outgoing_longwave_cubes, u_wind_850_cubes, u_wind_200_cubes,
                precip_cubes, runid, out_plot_dir)
            metrics.update(level3_metrics)
        '''
        if extreme_level4:
            extreme.precip_extremes_cs(precip_cubes, runid=runid, percentile=95, compute=True, plot=True)

        '''
        # write metrics to csv file for each suite_id
        with open(os.path.join(out_plot_dir, 'metrics.csv'), 'w') as fh:
            writer = csv.writer(fh)
            for metric in metrics.items():
                writer.writerow(metric)
'''
