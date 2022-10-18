import os, sys
import iris
import csv
# use Agg backend for running without X-server
import matplotlib as mpl
from mjo import data_paths
from mjo.src import diags_level1 as diags_level1
from mjo.src import diags_level2 as diags_level2
from mjo.src import diags_level3 as diags_level3
from mjo.src import diags_level4_prop as diags_level4_prop
from mjo.src import reg_lat_long_grid as reg_lat_long_grid

mpl.use('Agg')


def regrid_to_2p5degree_tropics(cube):
    """
    Regrids data to 2.5 x 2.5 degree grid and extract tropics

    :param cube: Initial cube
    :type cube:
    :return: Regridded cube
    :rtype:
    """
    # Extract region
    cube = cube.intersection(latitude=(-30, 30), longitude=(0, 360))
    # regrid often fails due to attributes mismatch!
    # A very annoying feature of Iris. So removing all attributes
    # cube.attributes = {}

    # regrid to regular lat-long 2.5 degree grid
    ref_cube = reg_lat_long_grid.create_cube(latitudes=(-30, 30),
                                             longitudes=(0, 360), spacing=2.5)

    # Copy the coordinate system over to the incoming cube
    cube.coord('latitude').coord_system = ref_cube.coord('latitude').coord_system
    cube.coord('longitude').coord_system = ref_cube.coord('longitude').coord_system

    if cube.coord('longitude').bounds is None:
        cube.coord('longitude').guess_bounds()
    if cube.coord('latitude').bounds is None:
        cube.coord('latitude').guess_bounds()
    return cube.regrid(ref_cube, iris.analysis.Linear())


def mjo_compute(control=None, expt=None, obs=None,
                level1=True,
                level2=True,
                level3=True,
                level4_prop=True):
    '''
    Compute the MJO diagnostics
    MJO computation is now part of SEAPy.

    It has 4 sections:
    1. Mean and variance plots
    2. Wavenumber-Frequency spectra
    3. Wheeler-Hendon plots and phase composites
    4. Propagation plots and statistics

    :param control: Baseline model
    :type control: Dictionary
    :param expt: Experiment model
    :type expt: Dictionary
    :param obs: Observations data
    :type obs: Dictionary
    :param level1: Compute level 1 diagnostics (Mean and variance plots)
    :type level1: Logical
    :param level2: Compute level 2 diagnostics (Wavenumber-Frequency spectra)
    :type level2: Logical
    :param level3: Compute level 3 diagnostics (Wheeler-Hendon plots and phase composites)
    :type level3: Logical
    :param level4_prop: Compute level 4 diagnostics (Propagation plots and statistics)
    :type level4_prop: Logical
    :return:
    :rtype:
    '''

    runs = [control, expt, obs]

    # Pick up only 'non-None' (valid) runs
    runs = [r for r in runs if r != None]

    for run in runs:
        metrics = {}
        runid = run['runid']
        label = run['label']
        # data
        data_root = os.path.join(run['data_retrieve_dir'], runid)

        # Out Plot Directory
        out_plot_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
        if not os.path.exists(out_plot_dir):
            os.makedirs(out_plot_dir)

        varnames = ['OLR', 'U850', 'U200', 'PRECIP', 'SST']

        for varname in varnames:
            combined_ppfile = os.path.join(data_root, runid + '_' + varname + '.pp.nc')
            # print(combined_ppfile)

            if varname == 'OLR':
                outgoing_longwave_cubes = iris.load_cube(combined_ppfile)
                outgoing_longwave_cubes.long_name = 'toa_outgoing_longwave_flux'
                outgoing_longwave_cubes = outgoing_longwave_cubes.intersection(latitude=(-30, 30),
                                                                               longitude=(0, 360))
                outgoing_longwave_cubes = regrid_to_2p5degree_tropics(outgoing_longwave_cubes)

            if varname == 'U850':
                u_wind_850_cubes = iris.load_cube(combined_ppfile)
                u_wind_850_cubes.long_name = 'x_wind_850hPa'
                u_wind_850_cubes = u_wind_850_cubes.intersection(latitude=(-30, 30),
                                                                 longitude=(0, 360))
                u_wind_850_cubes = regrid_to_2p5degree_tropics(u_wind_850_cubes)

            if varname == 'U200':
                u_wind_200_cubes = iris.load_cube(combined_ppfile)
                u_wind_200_cubes.long_name = 'x_wind_200hPa'
                u_wind_200_cubes = u_wind_200_cubes.intersection(latitude=(-30, 30),
                                                                 longitude=(0, 360))
                u_wind_200_cubes = regrid_to_2p5degree_tropics(u_wind_200_cubes)

            if varname == 'PRECIP':
                precip_cubes = iris.load_cube(combined_ppfile)
                precip_cubes.long_name = 'precipitation_flux'
                precip_cubes = precip_cubes.intersection(latitude=(-30, 30),
                                                         longitude=(0, 360))
                precip_cubes = regrid_to_2p5degree_tropics(precip_cubes)
                # for cube in precip_cubes:
                # precip_cubes.convert_units('kg m-2 day-1')
                if not runid == 'obs':
                    precip_cubes.convert_units('kg m-2 day-1')
            '''
            if varname == 'SST':
                sst_cubes = iris.load_cube(combined_ppfile)
                sst_cubes.long_name = 'surface_temperature'
                sst_cubes = sst_cubes.intersection(latitude=(-30, 30),
                                                         longitude=(0, 360))
                sst_cubes = regrid_to_2p5degree_tropics(sst_cubes)
            '''
        print(outgoing_longwave_cubes, u_wind_200_cubes, u_wind_850_cubes)

        # Level 1 diagnostics
        # Mean, variance, filtered variance, filt variance/total variance
        if level1:
            for cube in [precip_cubes, outgoing_longwave_cubes, u_wind_200_cubes,
                         u_wind_850_cubes]:
                diags_level1.diagnos_level1(cube, out_plot_dir, runid, label)

        # Level 2 diagnostics
        # WK raw sym/antisym spectra, background spectra, (sym, antisym)/background
        if level2:
            for cube in [precip_cubes, outgoing_longwave_cubes, u_wind_200_cubes,
                         u_wind_850_cubes]:
            #for cube in [precip_cubes]:
                print(cube)
                level2_metrics = diags_level2.diagnos_level2(cube, out_plot_dir, runid, label)
                metrics.update(level2_metrics)

        # Level 3 diagnostics
        # Real-time multivariate MJO Index (RMM) calculations, and Summer/Winter
        # composites
        if level3:
            level3_metrics = diags_level3.diagnos_level3(
                outgoing_longwave_cubes, u_wind_850_cubes, u_wind_200_cubes,
                precip_cubes, runid, label, out_plot_dir)
            metrics.update(level3_metrics)

        if level4_prop:
            diags_level4_prop.diagnos_level4_prop(outgoing_longwave_cubes, u_wind_850_cubes, u_wind_200_cubes,
                                                  precip_cubes, run, out_plot_dir=out_plot_dir)

        # write metrics to csv file for each suite_id
        with open(os.path.join(out_plot_dir, 'metrics.csv'), 'w') as fh:
            writer = csv.writer(fh)
            for metric in metrics.items():
                writer.writerow(metric)


