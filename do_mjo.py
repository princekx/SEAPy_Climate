"""MJO diagnostics driver.

Mirrors the layout used by :mod:`do_sea`: every run writes its plots and small
NetCDF files to ``<mjo_plots_dir>/<runid>_<dir_label>/<runid>_*``.
"""

import csv
import logging
import os
import warnings

import matplotlib
matplotlib.use('Agg', force=True)

import iris
from iris.warnings import IrisUserWarning, IrisVagueMetadataWarning

if hasattr(iris, 'FUTURE'):
    iris.FUTURE.date_microseconds = True
    iris.FUTURE.save_split_attrs = True

warnings.filterwarnings('ignore', category=IrisVagueMetadataWarning)
warnings.filterwarnings('ignore', category=IrisUserWarning)
warnings.filterwarnings('ignore', category=Warning, module=r'iris$')
warnings.filterwarnings('ignore', category=Warning, module=r'iris\..*')

from src import diags_mjo_mean
from src import diags_mjo_spectra
from src import diags_mjo_composites
from src import diags_level4_prop
from src import reg_lat_long_grid


logger = logging.getLogger(__name__)

MJO_VARIABLES = {
    'OLR': 'toa_outgoing_longwave_flux',
    'U850': 'x_wind_850hPa',
    'U200': 'x_wind_200hPa',
    'PRECIP': 'precipitation_flux',
}

MJO_VARNAMES = [
    'precipitation_flux',
    'toa_outgoing_longwave_flux',
    'x_wind_200hPa',
    'x_wind_850hPa',
]

SEASONS = ['summer', 'winter']


def _run_output_label(run):
    dir_label = run['dir_label'].replace(' ', '_')
    return f"{run['runid']}_{dir_label}"


def _run_plot_label(run):
    if run['runid'] == run.get('model'):
        return f"{run['runid']} \n {run['label']}"
    return f"{run['runid']} \n {run.get('model', '')} {run['label']}".replace('  ', ' ')


def _resolve_data_root(run):
    base_dir = run['data_retrieve_dir']
    dir_label = run['dir_label'].replace(' ', '_')
    runid = run['runid']

    candidates = [
        os.path.join(base_dir, f'{runid}_{dir_label}'),
        os.path.join(base_dir, f'{dir_label}_{runid}'),
        os.path.join(base_dir, dir_label),
        os.path.join(base_dir, runid),
        base_dir,
    ]

    for candidate in candidates:
        if os.path.isdir(candidate):
            logger.info('Using data directory %s for %s', candidate, runid)
            return candidate

    raise FileNotFoundError(
        f"Could not locate data directory for {runid}. Tried: {', '.join(candidates)}"
    )


def expected_plot_files(run, out_plot_dir, level1=True, level2=True, level3=True,
                        level4_prop=True, format='png'):
    """Return the plot paths each enabled MJO diagnostic level is expected to produce."""
    runid = run['runid']
    expected = []

    if level1:
        for varname in MJO_VARNAMES:
            for product in ['Mean', 'Var', 'VarRatio']:
                for season in SEASONS:
                    expected.append(
                        os.path.join(out_plot_dir, f'{runid}_{varname}_{product}_{season}.{format}')
                    )

    if level2:
        for varname in MJO_VARNAMES:
            expected.append(os.path.join(out_plot_dir, f'{runid}_{varname}_Ratio_Spec_Asym.{format}'))
            expected.append(os.path.join(out_plot_dir, f'{runid}_{varname}_Ratio_Spec_Sym.{format}'))
            for season in SEASONS:
                expected.append(
                    os.path.join(
                        out_plot_dir,
                        f'{runid}_{varname}_wavenum_freq_season_{season}.{format}',
                    )
                )

    if level3:
        expected.append(os.path.join(out_plot_dir, f'{runid}_WH04.{format}'))
        expected.append(os.path.join(out_plot_dir, f'{runid}_phase_freq_plot.{format}'))
        for varname in MJO_VARNAMES:
            expected.append(os.path.join(out_plot_dir, f'{runid}_{varname}_LeadLagCorr.{format}'))
            for season in SEASONS:
                expected.append(
                    os.path.join(
                        out_plot_dir,
                        f'{runid}_{varname}_{season}_phase_composites.{format}',
                    )
                )

    if level4_prop:
        for varname in MJO_VARNAMES:
            expected.append(os.path.join(out_plot_dir, f'{runid}_{varname}_prop_comp.{format}'))
            expected.append(os.path.join(out_plot_dir, f'{runid}_{varname}_filt_prop_comp.{format}'))

    return expected


def missing_plot_files(run, out_plot_dir, **kwargs):
    """Return the expected MJO plots that are not yet on disk for a run."""
    return [
        path
        for path in expected_plot_files(run, out_plot_dir, **kwargs)
        if not os.path.exists(path)
    ]


def regrid_to_2p5degree_tropics(cube):
    """Regrid a cube onto a regular 2.5 degree tropical grid."""
    cube = cube.intersection(latitude=(-30, 30), longitude=(0, 360))

    ref_cube = reg_lat_long_grid.create_cube(latitudes=(-30, 30), longitudes=(0, 360), spacing=2.5)

    cube.coord('latitude').coord_system = ref_cube.coord('latitude').coord_system
    cube.coord('longitude').coord_system = ref_cube.coord('longitude').coord_system

    if cube.coord('longitude').bounds is None:
        cube.coord('longitude').guess_bounds()
    if cube.coord('latitude').bounds is None:
        cube.coord('latitude').guess_bounds()

    return cube.regrid(ref_cube, iris.analysis.Linear())


def _normalise_time_coord(cube):
    """Rename an Iris time axis called valid_time so downstream code finds 'time'."""
    if cube.coords('time'):
        return cube

    for coord in cube.coords():
        if coord.name() == 'valid_time' or coord.standard_name == 'time':
            coord.rename('time')
            return cube

    raise iris.exceptions.CoordinateNotFoundError(
        f"No usable time coordinate on cube {cube.name()}"
    )


def _load_run_cubes(run, data_root):
    """Load, rename and regrid the cubes needed by the MJO diagnostics."""
    runid = run['runid']
    cubes = {}

    for varname, long_name in MJO_VARIABLES.items():
        file_path = os.path.join(data_root, f'{runid}_{varname}.pp.nc')
        if not os.path.exists(file_path):
            logger.warning('Missing %s file for %s: %s', varname, runid, file_path)
            continue

        try:
            cube = iris.load_cube(file_path)
        except Exception:
            logger.exception('Failed to load %s from %s', varname, file_path)
            continue

        try:
            cube = _normalise_time_coord(cube)
        except iris.exceptions.CoordinateNotFoundError:
            logger.warning('Skipping %s for %s: no time coordinate in %s', varname, runid, file_path)
            continue

        cube.long_name = long_name
        cube = regrid_to_2p5degree_tropics(cube)

        if varname == 'PRECIP' and runid != 'obs':
            try:
                cube.convert_units('kg m-2 day-1')
            except Exception:
                logger.exception('Failed to convert %s precipitation units', runid)
                continue

        cubes[varname] = cube

    return cubes


def mjo_compute(control=None, expt=None, obs=None,
                level1=True, level2=True, level3=True, level4_prop=True,
                mjo_plots_dir=None):
    """Run the MJO diagnostics for each requested dataset.

    Levels:
      1. Mean and variance maps
      2. Wavenumber-frequency spectra
      3. Wheeler-Hendon index and phase composites
      4. Propagation composites and statistics

    A failure in one level is logged and does not stop the remaining levels.
    """
    if not mjo_plots_dir:
        raise ValueError('mjo_plots_dir must be provided')

    runs = [run for run in [control, expt, obs] if run is not None]

    for run in runs:
        metrics = {}

        try:
            runid = run['runid']
            label = _run_plot_label(run)
            out_plot_dir = os.path.join(mjo_plots_dir, _run_output_label(run))
            os.makedirs(out_plot_dir, exist_ok=True)

            missing_plots = missing_plot_files(
                run,
                out_plot_dir,
                level1=level1,
                level2=level2,
                level3=level3,
                level4_prop=level4_prop,
            )
            if not missing_plots:
                logger.info(
                    'All expected MJO plots already exist for %s in %s; skipping run without opening data',
                    runid,
                    out_plot_dir,
                )
                continue

            logger.info(
                '%s expected MJO plots missing for %s (e.g. %s); running diagnostics',
                len(missing_plots),
                runid,
                os.path.basename(missing_plots[0]),
            )

            data_root = _resolve_data_root(run)
            cubes = _load_run_cubes(run, data_root)

            required = set(MJO_VARIABLES)
            if not required.issubset(cubes):
                logger.warning(
                    'Skipping MJO diagnostics for %s because these inputs are missing: %s',
                    runid,
                    ', '.join(sorted(required - set(cubes))),
                )
                continue

            olr = cubes['OLR']
            u850 = cubes['U850']
            u200 = cubes['U200']
            precip = cubes['PRECIP']

            if level1:
                for cube in [precip, olr, u200, u850]:
                    try:
                        diags_mjo_mean.diagnos_level1(cube, out_plot_dir, runid, label)
                    except Exception:
                        logger.exception('Level 1 diagnostics failed for %s %s', runid, cube.long_name)

            if level2:
                for cube in [precip, olr, u200, u850]:
                    try:
                        metrics.update(diags_mjo_spectra.diagnos_level2(cube, out_plot_dir, runid))
                    except Exception:
                        logger.exception('Level 2 diagnostics failed for %s %s', runid, cube.long_name)

            if level3 or level4_prop:
                olr, u850, u200, precip = diags_mjo_composites.check_common_time_period(
                    [olr, u850, u200, precip]
                )

            if level3:
                try:
                    metrics.update(
                        diags_mjo_composites.diagnos_level3(olr, u850, u200, precip, runid, label, out_plot_dir)
                    )
                except Exception:
                    logger.exception('Level 3 diagnostics failed for %s', runid)

            if level4_prop:
                try:
                    diags_level4_prop.diagnos_level4_prop(
                        olr, u850, u200, precip, run, out_plot_dir=out_plot_dir
                    )
                except Exception:
                    logger.exception('Level 4 propagation diagnostics failed for %s', runid)

            metrics_file = os.path.join(out_plot_dir, f'{runid}_metrics.csv')
            try:
                with open(metrics_file, 'w', newline='') as fh:
                    writer = csv.writer(fh)
                    for metric in metrics.items():
                        writer.writerow(metric)
                logger.info('Written metrics to %s', metrics_file)
            except Exception:
                logger.exception('Failed to write MJO metrics CSV for %s', runid)

        except Exception:
            logger.exception('Unhandled error while processing MJO run %s', run.get('runid', 'unknown'))
