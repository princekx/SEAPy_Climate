"""Equatorial-wave diagnostics.

This module computes seasonal variance, filtered wave composites and summary plots
for MJO, Kelvin and Rossby waves. All run-specific outputs are written into the run
folder implied by the YAML configuration, for example:

    /data/users/prince.xavier/SEAPy_climate/SEAPy_output/u-ea560_GC6.0_AMIP_C144/

The public entry point is :func:`eqwaves_compute`.
"""

import logging
import os

import cartopy.crs as ccrs
import iris
import iris.plot as iplt
import matplotlib
import matplotlib.pyplot as plt
from .plot_utils import save_figure
import numpy as np
import numpy.ma as ma
from iris.coord_categorisation import add_day_of_month, add_month_number, add_year
from matplotlib.colors import BoundaryNorm
from scipy.signal import find_peaks

from src import kf_filter

matplotlib.use('Agg')
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
    """Return the run label used for the directory name, matching the YAML config."""
    run = run or {}
    dir_label = str(run.get('dir_label', run.get('label', ''))).replace(' ', '_')
    return dir_label


def _resolve_wave_output_dir(run=None, out_plot_dir=None):
    """Resolve the output directory for a single run.

    The YAML-driven run output is expected to use the pattern
    ``<SEAPy_output>/<runid>_<dir_label>``. When ``out_plot_dir`` is provided, it is
    preferred and used as-is.
    """
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


def _output_path(run, filename, out_plot_dir=None):
    """Build a full output path inside the run-specific output directory."""
    out_dir = _resolve_wave_output_dir(run=run, out_plot_dir=out_plot_dir)
    _ensure_dir(out_dir)
    return os.path.join(out_dir, filename)


def _normalise_time_coord(cube):
    """Accept files whose time axis is named valid_time instead of time."""
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


def _add_time_coords(cube):
    """Ensure the time metadata needed for seasonal grouping is present."""
    cube = _normalise_time_coord(cube)
    if not cube.coords('year'):
        add_year(cube, 'time', name='year')
    if not cube.coords('month_number'):
        add_month_number(cube, 'time', name='month_number')
    if not cube.coords('day_of_month'):
        add_day_of_month(cube, 'time', name='day_of_month')
    return cube


def _contour_info(varname, wavename):
    """Return plotting levels for a given field and wave type."""
    info = {}

    if varname == 'precipitation_flux':
        info['shade_levels'] = np.arange(0, 44, 4)
        info['contour_levels'] = np.arange(0, 0.5, 0.05)
        if wavename == 'ALL':
            info['std_levels'] = np.arange(1, 40., 5)
        elif wavename in {'MJO', 'Kelvin'}:
            info['std_levels'] = np.arange(1, 4., 0.25)
            info['anom_levels'] = np.linspace(-4, 4., 11)
        else:
            info['std_levels'] = np.arange(1, 3., 0.25)

    elif varname == 'toa_outgoing_longwave_flux':
        info['shade_levels'] = np.arange(150, 240, 10)
        info['contour_levels'] = np.arange(150, 240, 10)
        if wavename == 'ALL':
            info['std_levels'] = np.arange(10, 70., 5)
        elif wavename in {'MJO', 'Kelvin'}:
            info['std_levels'] = np.arange(5, 16, 1)
            info['anom_levels'] = np.linspace(-15, 15, 11)
        else:
            info['std_levels'] = np.arange(5, 16, 1)

    elif varname in {'x_wind_850hPa', 'x_wind_200hPa'}:
        info['shade_levels'] = np.arange(0, 2.2, 0.2)
        info['contour_levels'] = np.arange(0, 2.2, 0.2)
        if wavename == 'ALL':
            info['std_levels'] = np.arange(1, 5, 0.25)
        elif wavename in {'MJO', 'Kelvin'}:
            info['std_levels'] = np.arange(1, 3., 0.25) if wavename == 'MJO' else np.arange(0.5, 1.5, 0.1)
            info['anom_levels'] = np.linspace(-2, 2., 11)
        else:
            info['std_levels'] = np.arange(0.5, 1.5, 0.1)

    elif varname == 'y_wind_850hPa':
        info['shade_levels'] = np.arange(0, 2.2, 0.2)
        info['contour_levels'] = np.arange(0, 2.2, 0.2)
        if wavename == 'ALL':
            info['std_levels'] = np.arange(1, 5, 0.25)
        else:
            info['std_levels'] = np.arange(0.1, 1., 0.1)
            info['anom_levels'] = np.linspace(-1, 1., 11)

    return info


def plot_wave_std_maps(wavename=None, varname=None, season=None, pltname='iris_test_plot.ps', run=None):
    """Plot a seasonal standard-deviation map for one wave type."""
    run = run or {}
    runid = run.get('runid', 'unknown')
    dir_label = _run_label(run)
    out_dir = _resolve_wave_output_dir(run=run)

    file_name = f'{runid}_{varname}_{wavename}_{season}_variance.nc'
    variance_file = os.path.join(out_dir, file_name)
    if os.path.exists(variance_file):
        cube = iris.load_cube(variance_file)
        cube.data = np.sqrt(cube.data)
    else:
        raise FileNotFoundError(f'Missing variance cube: {variance_file}')

    loni, lonf, lati, latf = (90, 140, -10, 20)
    cube = cube.intersection(latitude=(lati, latf), longitude=(loni, lonf))
    levels = _contour_info(cube.long_name, wavename).get('std_levels', np.arange(0.1, 2.0, 0.1))

    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    cmap = plt.get_cmap('YlOrRd')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cf = iplt.pcolormesh(cube, cmap=cmap, norm=norm)
    plt.colorbar(cf, orientation='horizontal')
    ax.gridlines(draw_labels=True, linewidth=0.25, alpha=0.5)
    ax.gridlines().xlabels_top = False
    ax.gridlines().ylabels_right = False
    plt.title(f'{dir_label} {cube.long_name} {wavename} amplitude {season}')
    ax.coastlines()
    save_figure(pltname, bbox_inches='tight', pad_inches=0)
    plt.close()


def plot_wave_variance_percent(wavename=None, varname=None, season=None, pltname='iris_test_plot.ps', run=None):
    """Plot the percentage of variance explained by a specific wave type."""
    run = run or {}
    runid = run.get('runid', 'unknown')
    dir_label = _run_label(run)
    out_dir = _resolve_wave_output_dir(run=run)

    variance_file = os.path.join(out_dir, f'{runid}_{varname}_{wavename}_{season}_variance.nc')
    if not os.path.exists(variance_file):
        raise FileNotFoundError(f'Missing variance cube: {variance_file}')
    cube = iris.load_cube(variance_file)

    all_variance_file = os.path.join(out_dir, f'{runid}_{varname}_ALL_{season}_variance.nc')
    if os.path.exists(all_variance_file):
        cube_all = iris.load_cube(all_variance_file)
        threshold = (cube_all.data.max() - cube_all.data.min()) * 0.01
        cube_all.data = ma.masked_less(cube_all.data, threshold)
        cube.data = 100.0 * cube.data / cube_all.data
    else:
        logger.warning('No ALL-wave variance file for %s (%s)', wavename, season)

    loni, lonf, lati, latf = (90, 140, -10, 20)
    cube = cube.intersection(latitude=(lati, latf), longitude=(loni, lonf))

    levels = [2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45]
    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    cmap = plt.get_cmap('ocean_r')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cf = iplt.pcolormesh(cube, cmap=cmap, norm=norm)
    plt.colorbar(cf, orientation='horizontal')
    ax.gridlines(draw_labels=True, linewidth=0.25, alpha=0.5)
    ax.gridlines().xlabels_top = False
    ax.gridlines().ylabels_right = False
    plt.title(f'{dir_label} {cube.long_name} {wavename} variance % {season}')
    ax.coastlines()
    save_figure(pltname, bbox_inches='tight', pad_inches=0)
    plt.close()


def plot_hov_composite(wavename=None, varname=None, run=None, pltname='iris_test_plot.ps'):
    """Plot a Hovmoller-style composite of wave anomalies along longitude."""
    run = run or {}
    runid = run.get('runid', 'unknown')
    dir_label = _run_label(run)
    out_dir = _resolve_wave_output_dir(run=run)
    composite_file = os.path.join(out_dir, f'{runid}_{varname}_{wavename}_refEquator_composite.nc')
    if not os.path.exists(composite_file):
        raise FileNotFoundError(f'Missing composite file: {composite_file}')

    cube = iris.load_cube(composite_file).intersection(longitude=(40, 180))
    if wavename == 'Kelvin':
        cube = cube.intersection(latitude=(-5, 5))
    else:
        cube = cube.intersection(latitude=(-10, 10))

    hov_cube = cube.collapsed('latitude', iris.analysis.MEAN)
    times = np.linspace(-len(hov_cube.coord('time').points) / 2, len(hov_cube.coord('time').points) / 2,
                        len(hov_cube.coord('time').points))
    lons = hov_cube.coord('longitude').points
    L, T = np.meshgrid(lons, times)

    levels = _contour_info(varname, wavename).get('anom_levels', np.linspace(-5, 5, 11))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cf = ax.contourf(L, T, hov_cube.data, levels=levels, cmap='RdBu', extend='both')
    fig.colorbar(cf, orientation='horizontal')
    plt.xlim(min(lons), max(lons))
    if wavename == 'Kelvin':
        plt.ylim(-10, 10)
    else:
        plt.ylim(-35, 35)
    plt.title(f'{dir_label} {varname} {wavename} composite')
    plt.xlabel('Longitude (degrees east)')
    plt.ylabel('Lag/lead (days)')
    plt.grid()
    save_figure(pltname, bbox_inches='tight')
    plt.close(fig)


def compute_wave_composite(cube, lat_range=(-10, 10), lon_range=(80, 100), wavename=None, varname=None, run=None):
    """Build a composite centred on wave events within a longitude-latitude box."""
    run = run or {}
    runid = run.get('runid', 'unknown')
    out_dir = _resolve_wave_output_dir(run=run)
    composite_file = os.path.join(out_dir, f'{runid}_{varname}_{wavename}_refEquator_composite.nc')

    cube = _normalise_time_coord(cube)
    ntime, nlat, nlon = cube.shape

    cube_ts = cube.intersection(latitude=lat_range, longitude=lon_range).collapsed(
        ['latitude', 'longitude'],
        iris.analysis.MEAN,
    )
    cube_ts_normal = (cube_ts - cube_ts.collapsed('time', iris.analysis.MEAN)) / cube_ts.collapsed(
        'time',
        iris.analysis.STD_DEV,
    )

    peaks, _ = find_peaks(cube_ts_normal.data, height=1.0)
    lag = 50
    comp = np.ndarray((len(peaks), 2 * lag + 1, nlat, nlon), dtype=float)

    for idx, peak in enumerate(peaks):
        if lag <= peak <= ntime - lag - 1:
            comp[idx] = cube[peak - lag: peak + lag + 1].data

    if len(peaks) == 0:
        logger.warning('No wave events found for %s on %s; writing empty composite.', wavename, runid)
        comp_cube = cube[: 2 * lag + 1].copy(data=np.zeros((2 * lag + 1, nlat, nlon), dtype=float))
    else:
        comp_cube = cube[: 2 * lag + 1].copy(data=np.mean(comp, axis=0))

    iris.save(comp_cube, composite_file)
    logger.info('Written %s', composite_file)
    return comp_cube


def compute_wave_variance_season(cube, run=None, wavename=None, metrics=None):
    """Compute seasonal variance and save the corresponding NetCDF files."""
    metrics = metrics or {}
    run = run or {}
    runid = run.get('runid', 'unknown')
    out_dir = _resolve_wave_output_dir(run=run)
    _ensure_dir(out_dir)

    cube = _normalise_time_coord(cube)
    cube = _add_time_coords(cube)
    varname = cube.long_name
    months = cube.coord('month_number').points

    ndjf_inds = [i for i in range(len(months)) if months[i] >= 11 or months[i] <= 2]
    jjas_inds = [i for i in range(len(months)) if 6 <= months[i] <= 9]

    ndjf_var = cube[ndjf_inds].collapsed('time', iris.analysis.VARIANCE)
    jjas_var = cube[jjas_inds].collapsed('time', iris.analysis.VARIANCE)

    ndjf_mean = ndjf_var.intersection(latitude=(-10, 12), longitude=(100, 130)).collapsed(
        ('latitude', 'longitude'),
        iris.analysis.MEAN,
    )
    jjas_mean = jjas_var.intersection(latitude=(-10, 12), longitude=(100, 130)).collapsed(
        ('latitude', 'longitude'),
        iris.analysis.MEAN,
    )

    metrics[f'Average NDJF {varname} {wavename} variance'] = round(float(ndjf_mean.data.item()), 2)
    metrics[f'Average JJAS {varname} {wavename} variance'] = round(float(jjas_mean.data.item()), 2)

    iris.save(ndjf_var, os.path.join(out_dir, f'{runid}_{varname}_{wavename}_ndjf_variance.nc'))
    iris.save(jjas_var, os.path.join(out_dir, f'{runid}_{varname}_{wavename}_jjas_variance.nc'))
    logger.info('Saved variance cubes for %s/%s', runid, wavename)
    return metrics


def eqwaves_compute(cube, out_plot_dir=None, run=None, compute=True, write_wave_data=True,
                    write_wave_data_dir=None, format='png'):
    """Run the equatorial-wave diagnostics for a single cube.

    The routine calculates the ALL, MJO, Kelvin and Rossby variance fields, writes
    them to run-specific NetCDF files, and saves the corresponding plots to the
    run-specific output folder.
    """
    run = run or {}
    metrics = {}
    cube = _normalise_time_coord(cube)
    cube = cube.intersection(latitude=(-30, 30), longitude=(0, 360))

    out_dir = _resolve_wave_output_dir(run=run, out_plot_dir=out_plot_dir)
    _ensure_dir(out_dir)
    if write_wave_data_dir is None:
        write_wave_data_dir = out_dir

    varname = cube.long_name or cube.name()
    runid = run.get('runid', 'unknown')
    dir_label = _run_label(run) or runid

    if not compute:
        logger.info('Wave computations are disabled; plotting only.')

    if compute:
        metrics.update(compute_wave_variance_season(cube.copy(), run=run, wavename='ALL', metrics=metrics))

    for season in ['ndjf', 'jjas']:
        path = os.path.join(out_dir, f'{runid}_{varname}_ALL_{season}_stddev.{format}')
        plot_wave_std_maps(wavename='ALL', varname=varname, season=season, pltname=path, run=run)

    for wavename in ['MJO', 'Kelvin', 'Rossby']:
        logger.info('Filtering %s...', wavename)
        if compute:
            ft = kf_filter.KFfilter(cube.data, spd=1)
            wave_cube = cube.copy()
            if wavename == 'MJO':
                wave_cube.data = ft.mjofilter()
            elif wavename == 'Kelvin':
                wave_cube.data = ft.kelvinfilter()
            else:
                wave_cube.data = ft.erfilter()

            try:
                compute_wave_composite(wave_cube, lat_range=(-10, 10), lon_range=(80, 100),
                                      wavename=wavename, varname=varname, run=run)
            except Exception:
                logger.exception('Failed to compute %s composite for %s', wavename, runid)

            metrics.update(compute_wave_variance_season(wave_cube, run=run, wavename=wavename, metrics=metrics))

            if write_wave_data:
                iris.save(wave_cube, os.path.join(write_wave_data_dir, f'{runid}_{varname}_{wavename}_filtered.nc'))

        for season in ['ndjf', 'jjas']:
            std_path = os.path.join(out_dir, f'{runid}_{varname}_{wavename}_{season}_stddev.{format}')
            var_path = os.path.join(out_dir, f'{runid}_{varname}_{wavename}_{season}_variance_percent.{format}')
            plot_wave_std_maps(wavename=wavename, varname=varname, season=season, pltname=std_path, run=run)
            plot_wave_variance_percent(wavename=wavename, varname=varname, season=season, pltname=var_path, run=run)

        try:
            hov_path = os.path.join(out_dir, f'{dir_label}_{varname}_{wavename}_composite_hovmoller.{format}')
            plot_hov_composite(wavename=wavename, varname=varname, run=run, pltname=hov_path)
        except Exception:
            logger.exception('Failed to plot %s hovmoller composite for %s', wavename, runid)

    logger.info('Wave diagnostics complete for %s', runid)
    return metrics


if __name__ == '__main__':
    obs = {
        'runid': 'obs',
        'dir_label': 'ERA5_GPM',
        'start_date': '1989/01/01',
        'end_date': '1998/12/01',
        'data_retrieve_dir': '/project/MJO_GCSS/hadgem3/data/obs/SEAPy_data',
    }
    print(f'Use eqwaves_compute(cube, run={obs}, out_plot_dir={_resolve_wave_output_dir(run=obs)})')

    run = obs
    data_root = os.path.join(run['data_retrieve_dir'], 'obs')
    out_plot_dir = '/project/MJO_GCSS/hadgem3/data/SEAPy_output'
    print(data_root)
    data_file = os.path.join(data_root, 'obs_PRECIP.pp.nc')
    cube = iris.load_cube(data_file)
    cube.long_name = 'precipitation_flux'
    #print(cube)

    wave_metrics = eqwaves_compute(cube, out_plot_dir, 'obs', compute=True,
                                   write_wave_data=True, write_wave_data_dir=data_root)
    # print(wave_metrics)
    #plot_hov_composite(wavename='MJO', varname='precipitation_flux', runid='obs')