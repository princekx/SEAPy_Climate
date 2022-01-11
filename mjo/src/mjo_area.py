#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
MJO assessment
"""
import argparse
import csv
import os
import os.path
import pickle
import tempfile

import iris
# use Agg backend for running without X-server
import matplotlib as mpl
mpl.use('Agg')

import mjo_utils as mu
import mjo_plots as mjp
import reg_lat_long_grid

import diags_level1
import diags_level2
import diags_level3
#import diags_level4


AUTOASSESS_PARSER = argparse.ArgumentParser(description='Description')
AUTOASSESS_PARSER.add_argument('--ref-suite-id', required=True, help='Reference suite ID.')
AUTOASSESS_PARSER.add_argument('--out-dir', required=True, help='Output directory')
AUTOASSESS_PARSER.add_argument('--data-dir', required=True, help='Directory containing model data.')
AUTOASSESS_PARSER.add_argument('--start-date', required=False, help='Start date for assessment.')
AUTOASSESS_PARSER.add_argument('--end-date', required=False, help='End date for assessment.')
AUTOASSESS_PARSER.add_argument('--obs-dir', required=False)
AUTOASSESS_PARSER.add_argument('--tmp-dir', required=True)
AUTOASSESS_PARSER.add_argument('--ancil-dir', required=False)
AUTOASSESS_PARSER.add_argument('--ncpu', default=1, help='Number of available processors.')


def parse_args(parser):
    """Parse command line arguments, and check all paths are absolute."""
    args = parser.parse_args()

    for arg, val in vars(args).items():
        if '-dir' in arg and not val is None and not os.path.isabs(val):
            msg = ' '.join(['Cli argument ', str(arg), ' not absolute path:', str(val)])
            raise NotAbsolutePath(msg)
    return args


def unpickle_cubes(path):
    """Load cube list from path."""
    with open(path, 'r') as fh:
        cubes = pickle.load(fh)
    return cubes


class NotAbsolutePath(Exception):
    pass


def main():
    AREA_NAME = 'MJO'.lower()

    args = parse_args(AUTOASSESS_PARSER)

    model_run_dirs = [
            os.path.join(args.data_dir, dir_) for dir_ in os.walk(args.data_dir).next()[1]
            ]

    cubes_paths = [
            os.path.join(dir_, AREA_NAME, 'cubes.pickle') for dir_ in model_run_dirs
            ]

    cubes = iris.cube.CubeList()
    for cubes_path in cubes_paths:
        cube_list = unpickle_cubes(cubes_path)
        cubes.extend(cube_list)

    # create output directories for single run assessments
    suite_ids = list(set(
        [cube.attributes['MODEL_RUN_ID_(added_by_AutoAssess)'] for cube in cubes]
        ))

    suite_ids.remove(args.ref_suite_id)
    suite_ids.append(args.ref_suite_id)
    AREA_OUT_DIR = os.path.join(args.out_dir, '_vs_'.join(suite_ids), AREA_NAME)

    for suite_id in suite_ids:
        out_dir = os.path.join(AREA_OUT_DIR, suite_id)
        os.makedirs(out_dir)

    # create tmp dir for large files
    tmp_dir = tempfile.mkdtemp(prefix='MJO_', dir=args.tmp_dir)


    ### Assessment
    # region (Tropics)
    lon = (0, 360)
    lat = (-30, 30)

    level1 = True
    level2 = True
    level3 = True
    level4 = True
    clean = True  # DEBUG option; delete temporary output

    # Extract region
    region_constraint = iris.Constraint(
            longitude=lambda cell: lon[0] <= cell <= lon[1],
            latitude=lambda cell: lat[0] <= cell <= lat[1]
            )
    cubes = cubes.extract(region_constraint)

    # regrid to regular lat-long 2.5 degree grid
    ref_cube = reg_lat_long_grid.create_cube(latitudes=(-30, 30),
            longitudes=(0, 360), spacing=2.5)

    regridded_cubes = iris.cube.CubeList()
    for cube in cubes:
        if cube.coord('longitude').bounds is None:
            cube.coord('longitude').guess_bounds()
        if cube.coord('latitude').bounds is None:
            cube.coord('latitude').guess_bounds()
        regridded_cubes.append(cube.regrid(ref_cube, iris.analysis.Linear()))
    cubes = regridded_cubes

    # assess single runs
    for suite_id in suite_ids:
        metrics = {}
        out_dir = os.path.join(AREA_OUT_DIR, suite_id)

        # TODO remove parentheses from attribute MODEL_RUN_ID_(added_by_AutoAssess)
        precip_cubes = cubes.extract(iris.Constraint('precipitation_flux'))
        for cube in precip_cubes:
            if cube.attributes['MODEL_RUN_ID_(added_by_AutoAssess)'] == suite_id:
                precip_cube = cube
                precip_cube.convert_units('kg m-2 day-1')

        outgoing_longwave_cubes = cubes.extract(iris.Constraint('toa_outgoing_longwave_flux'))
        for cube in outgoing_longwave_cubes:
            if cube.attributes['MODEL_RUN_ID_(added_by_AutoAssess)'] == suite_id:
                outgoing_longwave_cube = cube

        u_wind_200_cubes = cubes.extract(iris.Constraint('x_wind', pressure=200))
        for cube in u_wind_200_cubes:
            if cube.attributes['MODEL_RUN_ID_(added_by_AutoAssess)'] == suite_id:
                u_wind_200_cube = cube

        u_wind_850_cubes = cubes.extract(iris.Constraint('x_wind', pressure=850))
        for cube in u_wind_850_cubes:
            if cube.attributes['MODEL_RUN_ID_(added_by_AutoAssess)'] == suite_id:
                u_wind_850_cube = cube

        level1 = False
        # Level 1 diagnostics
        # Mean, variance, filtered variance, filt variance/total variance
        if level1:
            for cube in [precip_cube, outgoing_longwave_cube, u_wind_200_cube,
                    u_wind_850_cube]:

                diags_level1.diagnos_level1(cube, out_dir, suite_id, tmp_dir)

        level2 = True
        # Level 2 diagnostics
        # WK raw sym/antisym spectra, background spectra, (sym, antisym)/background
        if level2:
            for cube in [precip_cube, outgoing_longwave_cube, u_wind_200_cube,
                    u_wind_850_cube]:

                level2_metrics = diags_level2.diagnos_level2(cube, out_dir, suite_id)
                metrics.update(level2_metrics)

        level3 = False
        # Level 3 diagnostics
        # Real-time multivariate MJO Index (RMM) calculations, and Summer/Winter
        # composites
        if level3:
            level3_metrics = diags_level3.diagnos_level3(
                    outgoing_longwave_cube, u_wind_850_cube, u_wind_200_cube,
                    precip_cube, suite_id, out_dir
                    )
            metrics.update(level3_metrics)

        # write metrics to csv file for each suite_id
        with open(os.path.join(out_dir, 'metrics.csv'), 'w') as fh:
            writer = csv.writer(fh)
            for metric in metrics.items():
                writer.writerow(metric)

if __name__ == '__main__':
    main()
