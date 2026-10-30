"""Precipitation-extremes diagnostics.

These functions compute and plot percentile-based precipitation extremes for
cold-surge and MJO-driven periods. Outputs are written to the run-specific
SEAPy output directory so that a run such as:

    /data/users/prince.xavier/SEAPy_climate/SEAPy_output/u-ea560_GC6.0_AMIP_C144/

contains the relevant NetCDF and PNG products.
"""

import datetime
import logging
import os

import cartopy.crs as ccrs
import iris
import iris.plot as iplt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .plot_utils import save_figure
import numpy as np
import pandas as pd
from iris.coord_categorisation import add_day_of_month, add_month_number, add_year
from matplotlib.colors import BoundaryNorm

logger = logging.getLogger(__name__)

DEFAULT_OUTPUT_ROOT = '/data/users/prince.xavier/SEAPy_climate/SEAPy_output'


def _output_root():
    """Return the active output root, preferring the main workflow setting."""
    return os.environ.get('SEAPY_OUTPUT_DIR', DEFAULT_OUTPUT_ROOT)


def _ensure_dir(path):
    """Create a directory and return it."""
    os.makedirs(path, exist_ok=True)
    return path


def _run_label(run):
    """Return the directory label used by the run config."""
    if run is None:
        return ''
    return str(run.get('dir_label', run.get('label', ''))).replace(' ', '_')


def _resolve_output_dir(run=None, out_plot_dir=None):
    """Resolve the directory that should hold this run's precipitation plots."""
    if out_plot_dir:
        return out_plot_dir

    output_root = _output_root()
    if run is None:
        return output_root

    runid = str(run.get('runid', '')).strip()
    dir_label = _run_label(run)
    if runid and dir_label:
        return os.path.join(output_root, f'{runid}_{dir_label}')
    if dir_label:
        return os.path.join(output_root, dir_label)
    if runid:
        return os.path.join(output_root, runid)
    return output_root


def _normalise_time_coord(cube):
    """Accept a time axis stored as valid_time as well as time."""
    try:
        cube.coord('time')
        return cube
    except iris.exceptions.CoordinateNotFoundError:
        for coord in cube.coords():
            if coord.name() in {'time', 'valid_time'} or coord.standard_name == 'time':
                if coord.name() != 'time':
                    coord.rename('time')
                return cube
        raise


def prepare_calendar(cube):
    """Add year, month and day coordinates needed for date selection."""
    cube = _normalise_time_coord(cube)
    if not cube.coords('year'):
        add_year(cube, 'time', name='year')
    if not cube.coords('month_number'):
        add_month_number(cube, 'time', name='month_number')
    if not cube.coords('day_of_month'):
        add_day_of_month(cube, 'time', name='day_of_month')
    return cube


def return_indices_of_a(a, b):
    """Return positions in a where values match any value in b."""
    b_set = set(b)
    return [i for i, v in enumerate(a) if v in b_set]


def compute_percentiles(cube, run=None, cstype='NDJF', percent=(95,), cs_index_out_dir=None):
    """Compute percentile maps over a specified set of dates."""
    run = run or {}
    runid = run.get('runid', 'unknown')
    cs_index_out_dir = cs_index_out_dir or _resolve_output_dir(run=run)
    _ensure_dir(cs_index_out_dir)

    cs_dates_file = os.path.join(cs_index_out_dir, f'{runid}_{cstype}_dates.txt')
    if not os.path.exists(cs_dates_file):
        logger.warning('Date file not found: %s', cs_dates_file)
        return None

    df = pd.read_csv(cs_dates_file, header=None)
    df.columns = ['dates']
    selected_dates = [int(value) for value in df['dates']]

    cube = _normalise_time_coord(cube)
    times = cube.coord('time')
    cube_dates = [int(time.strftime('%Y%m%d')) for time in times.units.num2date(times.points)]
    indices = return_indices_of_a(cube_dates, selected_dates)

    if not indices:
        logger.warning('No matching dates found for %s in %s', cstype, runid)
        return None

    return cube[indices].collapsed('time', iris.analysis.PERCENTILE, percent=percent)


def precip_extremes_cs(precip_cubes, run=None, percentile=95,
                       compute=True, plot=True, frac_plot=True, out_plot_dir=None):
    """Compute, save and plot percentile-based precipitation extremes for CS seasons."""
    run = run or {}
    runid = run.get('runid', 'unknown')
    dir_label = _run_label(run)
    output_dir = _resolve_output_dir(run=run, out_plot_dir=out_plot_dir)
    _ensure_dir(output_dir)

    if compute:
        precip_cubes = _normalise_time_coord(precip_cubes)
        precip_cubes = precip_cubes.intersection(longitude=(95, 130), latitude=(-10, 25))
        precip_cubes = prepare_calendar(precip_cubes)

        for cst in ['NDJF', 'CS', 'CES']:
            precip_perc = compute_percentiles(
                precip_cubes,
                run=run,
                cstype=cst,
                percent=[percentile],
                cs_index_out_dir=output_dir,
            )
            if precip_perc is None:
                continue
            outfilename = os.path.join(output_dir, f'{runid}_{cst}_precip_{percentile}_percentile.nc')
            iris.save(precip_perc, outfilename)
            logger.info('Written %s', outfilename)

    if plot:
        for cst in ['NDJF', 'CS', 'CES']:
            filename = os.path.join(output_dir, f'{runid}_{cst}_precip_{percentile}_percentile.nc')
            if not os.path.exists(filename):
                continue
            precip_perc = iris.load_cube(filename)
            plt.figure()
            ax = plt.axes(projection=ccrs.PlateCarree())
            figname = os.path.join(output_dir, f'{runid}_{cst}_precip_{percentile}_percentile.png')

            levels = np.arange(10, 80, 10)
            cmap = plt.get_cmap('GnBu')
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            cf = iplt.pcolormesh(precip_perc, cmap=cmap, norm=norm)
            plt.colorbar(cf, orientation='vertical')
            plt.ylim([-10, 25])
            plt.xlim([95, 130])
            gl = ax.gridlines(draw_labels=True, alpha=0.5)
            gl.xlabels_top = False
            gl.ylabels_right = False
            plt.title(f'{dir_label} {cst} {percentile} percentile PRECIP')
            ax.coastlines()
            save_figure(figname, bbox_inches='tight', pad_inches=0)
            plt.close()
            logger.info('Plotted %s', figname)

    if frac_plot:
        seasonal_file = os.path.join(output_dir, f'{runid}_NDJF_precip_{percentile}_percentile.nc')
        if not os.path.exists(seasonal_file):
            logger.warning('Seasonal precipitation percentile file missing: %s', seasonal_file)
            return

        seasonal_precip_perc = iris.load_cube(seasonal_file)
        for cst in ['CS', 'CES']:
            filename = os.path.join(output_dir, f'{runid}_{cst}_precip_{percentile}_percentile.nc')
            if not os.path.exists(filename):
                continue
            precip_perc = iris.load_cube(filename)
            plt.figure()
            ax = plt.axes(projection=ccrs.PlateCarree())
            figname = os.path.join(output_dir, f'{runid}_{cst}_precip_{percentile}_percentile_fraction.png')
            levels = np.arange(-100, 120, 20)
            cmap = plt.get_cmap('RdBu')
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            data = 100.0 * (precip_perc - seasonal_precip_perc) / seasonal_precip_perc
            cf = iplt.pcolormesh(data, cmap=cmap, norm=norm)
            plt.colorbar(cf, orientation='vertical')
            plt.ylim([-10, 25])
            plt.xlim([95, 130])
            gl = ax.gridlines(draw_labels=True, alpha=0.5)
            gl.xlabels_top = False
            gl.ylabels_right = False
            plt.title(f'{dir_label} percent change in {percentile} percentile PRECIP due to {cst}')
            ax.coastlines()
            save_figure(figname, bbox_inches='tight', pad_inches=0)
            plt.close()
            logger.info('Plotted %s', figname)


def mjo_phase_indices(cube, run=None, season='NDJF', phases=None, mjo_index_out_dir=None):
    """Return cube indices matching the requested MJO phases within the season."""
    if phases is None:
        phases = [1]

    run = run or {}
    runid = run.get('runid', 'unknown')
    mjo_index_out_dir = mjo_index_out_dir or _resolve_output_dir(run=run)
    mjo_dates_file = os.path.join(mjo_index_out_dir, f'RMMs_{runid}.txt')
    if not os.path.exists(mjo_dates_file):
        raise FileNotFoundError(f'Missing MJO phase file: {mjo_dates_file}')

    df = pd.read_csv(mjo_dates_file, sep=' ', names=['year', 'month', 'day', 'rmm1', 'rmm2', 'phase', 'amp'])
    df['date'] = [datetime.datetime(y, m, d) for y, m, d in zip(df.year, df.month, df.day)]

    cube = _normalise_time_coord(cube)
    cube_dates = [int(time.strftime('%Y%m%d')) for time in cube.coord('time').units.num2date(cube.coord('time').points)]

    if season == 'NDJF':
        ndjf_df = df.loc[(df['month'].isin([11, 12, 1, 2]))]
        phase_df = ndjf_df.loc[(ndjf_df['amp'] >= 1) & (ndjf_df['phase'].isin(phases))]
        ndjf_dates = [int(date.strftime('%Y%m%d')) for date in ndjf_df['date']]
        phase_dates = [int(date.strftime('%Y%m%d')) for date in phase_df['date']]
        return return_indices_of_a(cube_dates, ndjf_dates), return_indices_of_a(cube_dates, phase_dates)

    raise ValueError(f'Unsupported season: {season}')


def precip_extremes_mjo(cube, run=None, percentile=95,
                        compute=True, plot=True, frac_plot=True, mjo_index_out_dir=None):
    """Compute percentile maps for MJO-phase precipitation extremes."""
    run = run or {}
    runid = run.get('runid', 'unknown')
    dir_label = _run_label(run)
    mjo_index_out_dir = mjo_index_out_dir or _resolve_output_dir(run=run)
    _ensure_dir(mjo_index_out_dir)

    cube = _normalise_time_coord(cube)
    cube = cube.intersection(longitude=(95, 130), latitude=(-10, 25))

    if compute:
        ndjf_indices, p234_indices = mjo_phase_indices(cube, phases=[2, 3, 4], run=run, season='NDJF',
                                                     mjo_index_out_dir=mjo_index_out_dir)
        seas_perc = cube[ndjf_indices].collapsed('time', iris.analysis.PERCENTILE, percent=[percentile])
        outfilename = os.path.join(mjo_index_out_dir, f'{runid}_NDJF_precip_{percentile}_percentile.nc')
        iris.save(seas_perc, outfilename)
        logger.info('Written %s', outfilename)

        mjo_perc_234 = cube[p234_indices].collapsed('time', iris.analysis.PERCENTILE, percent=[percentile])
        outfilename = os.path.join(mjo_index_out_dir, f'{runid}_MJO_234_precip_{percentile}_percentile.nc')
        iris.save(mjo_perc_234, outfilename)
        logger.info('Written %s', outfilename)

        ndjf_indices, p678_indices = mjo_phase_indices(cube, phases=[6, 7, 8], run=run, season='NDJF',
                                                     mjo_index_out_dir=mjo_index_out_dir)
        mjo_perc_678 = cube[p678_indices].collapsed('time', iris.analysis.PERCENTILE, percent=[percentile])
        outfilename = os.path.join(mjo_index_out_dir, f'{runid}_MJO_678_precip_{percentile}_percentile.nc')
        iris.save(mjo_perc_678, outfilename)
        logger.info('Written %s', outfilename)

    if plot:
        for mjo_phase in ['NDJF', 'MJO_234', 'MJO_678']:
            outfilename = os.path.join(mjo_index_out_dir, f'{runid}_{mjo_phase}_precip_{percentile}_percentile.nc')
            if not os.path.exists(outfilename):
                continue
            plot_cube = iris.load_cube(outfilename)
            plt.figure()
            ax = plt.axes(projection=ccrs.PlateCarree())
            figname = os.path.join(mjo_index_out_dir, f'{runid}_{mjo_phase}_precip_{percentile}_percentile.png')
            levels = np.arange(10, 80, 10)
            cmap = plt.get_cmap('GnBu')
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            cf = iplt.pcolormesh(plot_cube, cmap=cmap, norm=norm)
            plt.colorbar(cf, orientation='vertical')
            plt.ylim([-10, 25])
            plt.xlim([95, 130])
            gl = ax.gridlines(draw_labels=True, alpha=0.5)
            gl.xlabels_top = False
            gl.ylabels_right = False
            plt.title(f'{dir_label} {mjo_phase} {percentile} percentile PRECIP')
            ax.coastlines()
            save_figure(figname, bbox_inches='tight', pad_inches=0)
            plt.close()
            logger.info('Plotted %s', figname)

    if frac_plot:
        seasonal_file = os.path.join(mjo_index_out_dir, f'{runid}_NDJF_precip_{percentile}_percentile.nc')
        if not os.path.exists(seasonal_file):
            logger.warning('Seasonal precipitation file missing for MJO fraction plot: %s', seasonal_file)
            return
        seasonal_precip_perc = iris.load_cube(seasonal_file)
        for mjo_phase in ['MJO_234', 'MJO_678']:
            outfilename = os.path.join(mjo_index_out_dir, f'{runid}_{mjo_phase}_precip_{percentile}_percentile.nc')
            if not os.path.exists(outfilename):
                continue
            plot_cube = iris.load_cube(outfilename)
            plt.figure()
            ax = plt.axes(projection=ccrs.PlateCarree())
            figname = os.path.join(mjo_index_out_dir, f'{runid}_{mjo_phase}_precip_{percentile}_percentile_fraction.png')
            levels = np.arange(-100, 120, 20)
            cmap = plt.get_cmap('RdBu')
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            data = 100.0 * (plot_cube - seasonal_precip_perc) / seasonal_precip_perc
            cf = iplt.pcolormesh(data, cmap=cmap, norm=norm)
            plt.colorbar(cf, orientation='vertical')
            plt.ylim([-10, 25])
            plt.xlim([95, 130])
            gl = ax.gridlines(draw_labels=True, alpha=0.5)
            gl.xlabels_top = False
            gl.ylabels_right = False
            plt.title(f'{dir_label} percent change in {percentile} percentile PRECIP due to {mjo_phase}')
            ax.coastlines()
            save_figure(figname, bbox_inches='tight', pad_inches=0)
            plt.close()
            logger.info('Plotted %s', figname)


if __name__ == '__main__':
    sample_run = {
        'runid': 'obs',
        'dir_label': 'ERA5_GPM',
    }
    print(_resolve_output_dir(run=sample_run))
