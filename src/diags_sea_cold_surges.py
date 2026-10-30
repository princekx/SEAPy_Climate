import csv
import logging
import os
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import matplotlib
matplotlib.use('Agg', force=True)

import cartopy.crs as ccrs
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from .plot_utils import save_figure
import numpy as np
from iris.coord_categorisation import add_day_of_month, add_month_number, add_year
from iris.warnings import IrisUserWarning, IrisVagueMetadataWarning
from matplotlib.colors import BoundaryNorm

warnings.filterwarnings('ignore', category=IrisVagueMetadataWarning)
warnings.filterwarnings('ignore', category=IrisUserWarning)
warnings.filterwarnings('ignore', category=Warning, module=r'iris$')
warnings.filterwarnings('ignore', category=Warning, module=r'iris\..*')
warnings.filterwarnings(
    'ignore',
    message='.*legacy-style attribute handling.*',
    category=Warning,
)

if hasattr(iris, 'FUTURE'):
    iris.FUTURE.date_microseconds = True
    iris.FUTURE.save_split_attrs = True

logger = logging.getLogger(__name__)

DEFAULT_OUTPUT_ROOT = '/data/users/prince.xavier/SEAPy_climate/SEAPy_output'


def _output_root():
    """Return the active root output directory, preferring the workflow setting."""
    return os.environ.get('SEAPY_OUTPUT_DIR', DEFAULT_OUTPUT_ROOT)


def _resolve_output_dir(run=None, out_plot_dir=None):
    """Resolve the run-specific output directory for cold-surge products."""
    if out_plot_dir:
        return out_plot_dir

    run = run or {}
    runid = str(run.get('runid', '')).strip()
    dir_label = str(run.get('dir_label', run.get('label', ''))).replace(' ', '_')
    output_root = _output_root()
    if runid and dir_label:
        return os.path.join(output_root, f'{runid}_{dir_label}')
    if dir_label:
        return os.path.join(output_root, dir_label)
    if runid:
        return os.path.join(output_root, runid)
    return output_root


def _ensure_dir(path):
    """Create a directory tree if needed and return it."""
    os.makedirs(path, exist_ok=True)
    return path


def _write_date_file(file_path, dates):
    """Write a one-date-per-line text file using YYYYMMDD integers."""
    values = [str(int(date)) for date in dates]
    with open(file_path, 'w', encoding='utf-8') as handle:
        if values:
            handle.write('\n'.join(values) + '\n')


def _read_date_file(file_path):
    """Read a plain-text YYYYMMDD date list back into integers."""
    if not os.path.exists(file_path):
        return []
    with open(file_path, 'r', encoding='utf-8') as handle:
        return [int(line.strip()) for line in handle if line.strip()]


def _normalise_time_coord(cube):
    """Accept Iris files whose time axis is called valid_time instead of time."""
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
    """Ensure that year, month and day coordinates are present for a time axis."""
    cube = _normalise_time_coord(cube)
    if not cube.coords('year'):
        add_year(cube, 'time', name='year')
    if not cube.coords('month_number'):
        add_month_number(cube, 'time', name='month_number')
    if not cube.coords('day_of_month'):
        add_day_of_month(cube, 'time', name='day_of_month')
    return cube


def _cube_date_keys(cube):
    """Convert the time coordinate to integer YYYYMMDD keys for date matching."""
    cube = _add_time_coords(cube)
    years = cube.coord('year').points
    months = cube.coord('month_number').points
    days = cube.coord('day_of_month').points
    return np.asarray([int(years[i] * 10000 + months[i] * 100 + days[i]) for i in range(len(years))], dtype=int)


def _load_cube_if_exists(file_path):
    """Load a NetCDF cube if the file exists, otherwise return None."""
    if not os.path.exists(file_path):
        return None
    try:
        return iris.load_cube(file_path)
    except Exception:
        logger.exception('Failed to load cube from %s', file_path)
        return None


def _box_mean(cube, box_type):
    """Average a field over the specified latitude-longitude box."""
    subset = cube.intersection(latitude=(box_type[2], box_type[3]), longitude=(box_type[0], box_type[1]))
    return subset.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)


def prepare_data(u850, v850, box_type):
    """Compute box-mean wind and wind-speed time series for a selected region."""
    for cube in (u850, v850):
        _add_time_coords(cube)

    # Use the wind-speed magnitude directly rather than recalculating it repeatedly.
    speed_cube = u850.copy(data=np.hypot(u850.data, v850.data))
    u850_ba = _box_mean(u850, box_type)
    v850_ba = _box_mean(v850, box_type)
    speed_ba = _box_mean(speed_cube, box_type)
    return u850_ba, v850_ba, speed_ba.data


def consecutive_cs_days(cs_points):
    """Return lengths of contiguous cold-surge runs longer than 2 days."""
    values = np.asarray(cs_points, dtype=int)
    lengths = []
    current = 0
    for value in values:
        if value == 1:
            current += 1
        elif current > 2:
            lengths.append(current)
            current = 0
        else:
            current = 0

    if current > 2:
        lengths.append(current)

    return np.asarray(lengths, dtype=float)


def print_dict(metrics, run, cs_index_out_dir):
    """Write a metrics dictionary to a CSV file for a given run."""
    runid = run.get('runid', 'unknown')
    _ensure_dir(cs_index_out_dir)
    metrics_csv_file = os.path.join(cs_index_out_dir, f'{runid}_metrics.csv')
    try:
        with open(metrics_csv_file, 'w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerows(metrics.items())
    except Exception:
        logger.exception('Failed to write metrics CSV for %s', runid)
        raise

    logger.info('Written %s', metrics_csv_file)


def get_yyyymmdd(cube):
    """Return the time axis as integer YYYYMMDD values."""
    cube = _add_time_coords(cube)
    years = cube.coord('year').points
    months = cube.coord('month_number').points
    days = cube.coord('day_of_month').points
    return [int(years[i] * 10000 + months[i] * 100 + days[i]) for i in range(len(years))]


def find_array_match_inds(long_array, short_array):
    """Find insertion indices for matching values in a sorted array."""
    return np.asarray(long_array).searchsorted(np.asarray(short_array))


def return_indices_of_a(a, b):
    """Return positions in a where any member of b matches."""
    a_arr = np.asarray(a)
    b_set = set(int(v) for v in b)
    return np.nonzero(np.isin(a_arr, list(b_set)))[0].tolist()


def _mean_or_nan(values):
    """Mean of a numeric array, returning NaN when there are no valid values."""
    return float(np.mean(values)) if values.size else np.nan


def _count_per_year(values, years):
    """Average number of detected events per season, scaled by the total years in the series."""
    if years.size == 0:
        return np.nan
    season_years = int(years[-1] - years[0] + 1)
    return float(len(values)) / float(season_years) if season_years else np.nan


def cold_surge_stats(u850, v850, run=None, out_plot_dir=None):
    """Compute cold-surge statistics and write date-index files for each event class.

    The routine identifies NDJF dates, cold surges, cross-equatorial surges,
    meridional surges and easterly surges using simple box-based thresholds.
    """
    try:
        run = run or {}
        runid = run.get('runid', 'unknown')
        cs_index_out_dir = _ensure_dir(_resolve_output_dir(run=run, out_plot_dir=out_plot_dir))

        chang_box = [107, 115, 5, 10]
        chang_threshold = 9.0
        cs_threshold = 0.75

        u850_ba, v850_ba, speed_ba_ts = prepare_data(u850, v850, chang_box)
        years = _add_time_coords(u850_ba).coord('year').points
        months = u850_ba.coord('month_number').points
        date_keys = _cube_date_keys(u850_ba)

        ndjf_mask = (months >= 11) | (months <= 2)
        _write_date_file(os.path.join(cs_index_out_dir, f'{runid}_NDJF_dates.txt'), date_keys[ndjf_mask])

        csinds = np.where((speed_ba_ts >= chang_threshold) & (u850_ba.data < 0.0) & (v850_ba.data < 0.0))[0]
        cs_points = np.zeros(speed_ba_ts.shape, dtype=int)
        cs_points[csinds] = 1
        _write_date_file(os.path.join(cs_index_out_dir, f'{runid}_CS_dates.txt'), date_keys[csinds])

        hattori_box = [105, 115, -5, 5]
        hattori_threshold = -2.0
        u850_ha, v850_ha, _ = prepare_data(u850, v850, hattori_box)
        cesinds = np.where(
            (speed_ba_ts >= cs_threshold)
            & (u850_ba.data < 0.0)
            & (v850_ba.data < 0.0)
            & (v850_ha.data <= hattori_threshold)
        )[0]
        ces_points = np.zeros(speed_ba_ts.shape, dtype=int)
        ces_points[cesinds] = 1
        _write_date_file(os.path.join(cs_index_out_dir, f'{runid}_CES_dates.txt'), date_keys[cesinds])

        ms_box = [110, 117.5, 14, 16]
        ms_es_threshold = -8.0
        _, v850_ms, _ = prepare_data(u850, v850, ms_box)
        msinds = np.where(v850_ms.data < ms_es_threshold)[0]
        ms_points = np.zeros(v850_ms.shape, dtype=int)
        ms_points[msinds] = 1
        _write_date_file(os.path.join(cs_index_out_dir, f'{runid}_MS_dates.txt'), date_keys[msinds])

        es_box = [119, 121, 7.5, 15]
        u850_es, _, _ = prepare_data(u850, v850, es_box)
        esinds = np.where(u850_es.data < ms_es_threshold)[0]
        es_points = np.zeros(u850_es.shape, dtype=int)
        es_points[esinds] = 1
        _write_date_file(os.path.join(cs_index_out_dir, f'{runid}_ES_dates.txt'), date_keys[esinds])

        cs_stats = consecutive_cs_days(cs_points)
        ces_stats = consecutive_cs_days(ces_points)
        ms_stats = consecutive_cs_days(ms_points)
        es_stats = consecutive_cs_days(es_points)

        metrics = {
            'Average length of Cold Surges (CS)': round(_mean_or_nan(cs_stats), 2),
            'Average length of Cross-Equatorial Surges (CES)': round(_mean_or_nan(ces_stats), 2),
            'Average length of Meridional Surges (MS)': round(_mean_or_nan(ms_stats), 2),
            'Average length of Easterly Surges (ES)': round(_mean_or_nan(es_stats), 2),
            'Average number of CS in a season': round(_count_per_year(cs_stats, years), 2),
            'Average number of CES in a season': round(_count_per_year(ces_stats, years), 2),
            'Average number of MS in a season': round(_count_per_year(ms_stats, years), 2),
            'Average number of ES in a season': round(_count_per_year(es_stats, years), 2),
        }

        print(f'Output written to {cs_index_out_dir}')
        print_dict(metrics, run, cs_index_out_dir)
        return metrics
    except Exception:
        logger.exception('Cold surge stats failed for %s', run.get('runid', 'unknown'))
        raise


def _composite_one_variable(cube, cstype, runid, cs_index_out_dir):
    """Compute the composites for a single variable, returning (outfile, cube) pairs."""
    subset = cube.intersection(latitude=(-10, 25), longitude=(95, 130))
    all_dates = _cube_date_keys(subset)
    var_name = subset.long_name or subset.name()
    results = []

    for cst in cstype:
        logger.info('Computing composite of %s for %s', cst, var_name)
        cs_dates_file = os.path.join(cs_index_out_dir, f'{runid}_{cst}_dates.txt')
        if not os.path.exists(cs_dates_file):
            continue

        cs_dates = np.asarray(_read_date_file(cs_dates_file), dtype=int)
        if cs_dates.size == 0:
            continue

        inds = return_indices_of_a(all_dates, cs_dates)
        if not inds:
            logger.warning('No matching dates found for %s in %s; skipping composite.', cst, runid)
            continue

        comp = subset[inds].collapsed('time', iris.analysis.MEAN)
        _ = comp.data  # realise in the worker thread so only the write stays serial
        outfile = os.path.join(cs_index_out_dir, f'{runid}_{cst}_{var_name}_composite.nc')
        results.append((outfile, comp))

    return results


def cold_surge_composites(vars, cstype=None, run=None, out_plot_dir=None, max_workers=None):
    """Build time-mean composites for each identified cold-surge category.

    Each input cube is restricted to the Southeast Asia domain and then averaged
    over the dates listed in the corresponding *_dates.txt files. Variables are
    composited in parallel threads; the NetCDF writes stay serial because the
    underlying HDF5 library is not thread-safe.
    """
    try:
        run = run or {}
        runid = run.get('runid', 'unknown')
        cstype = cstype or ['NDJF', 'CS', 'CES', 'MS', 'ES']
        cs_index_out_dir = _ensure_dir(_resolve_output_dir(run=run, out_plot_dir=out_plot_dir))

        cubes = list(vars)
        if not cubes:
            return

        workers = max_workers or min(len(cubes), os.cpu_count() or 1)
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(_composite_one_variable, cube, cstype, runid, cs_index_out_dir): cube
                for cube in cubes
            }
            for future in as_completed(futures):
                for outfile, comp in future.result():
                    iris.save(comp, outfile)
                    logger.info('Written %s', outfile)
    except Exception:
        logger.exception('Cold surge composites failed for %s', run.get('runid', 'unknown'))
        raise


def _load_composite_file(cs_index_out_dir, runid, cst, var_name):
    """Helper for loading a previously written composite on demand."""
    return _load_cube_if_exists(os.path.join(cs_index_out_dir, f'{runid}_{cst}_{var_name}_composite.nc'))


def plot_cold_surge_composites(cstype=None, run=None, out_plot_dir=None, format='png'):
    """Plot precipitation and 850 hPa winds for composite cold-surge anomalies.

    The NDJF background state is used as the reference for anomaly composites,
    then each event category is displayed with a common Southeast Asia map layout.
    """
    try:
        run = run or {}
        runid = run.get('runid', 'unknown')
        cstype = cstype or ['CS', 'CES', 'MS', 'ES']
        dir_label = str(run.get('dir_label', run.get('label', ''))).replace(' ', '_')
        cs_index_out_dir = _ensure_dir(_resolve_output_dir(run=run, out_plot_dir=out_plot_dir))

        ndjf_precip = _load_composite_file(cs_index_out_dir, runid, 'NDJF', 'precipitation_flux')
        ndjf_sst = _load_composite_file(cs_index_out_dir, runid, 'NDJF', 'surface_temperature')
        ndjf_u850 = _load_composite_file(cs_index_out_dir, runid, 'NDJF', 'x_wind_850hPa')
        ndjf_v850 = _load_composite_file(cs_index_out_dir, runid, 'NDJF', 'y_wind_850hPa')

        if ndjf_sst.data.max() > 100:
            ndjf_sst.data -= 273.15

        if ndjf_precip is None or ndjf_u850 is None or ndjf_v850 is None:
            logger.warning('Missing NDJF composites for %s; skipping plot generation', runid)
            return

        for cst in cstype:
            if cst == 'NDJF':
                cst_precip_mean = ndjf_precip
                cst_sst_mean = ndjf_sst
                cst_u850_mean = ndjf_u850
                cst_v850_mean = ndjf_v850
            else:
                cst_precip_mean = _load_composite_file(cs_index_out_dir, runid, cst, 'precipitation_flux')
                cst_sst_mean = _load_composite_file(cs_index_out_dir, runid, cst, 'surface_temperature')
                cst_u850_mean = _load_composite_file(cs_index_out_dir, runid, cst, 'x_wind_850hPa')
                cst_v850_mean = _load_composite_file(cs_index_out_dir, runid, cst, 'y_wind_850hPa')

                if cst_precip_mean is None or cst_sst_mean is None or cst_u850_mean is None or cst_v850_mean is None:
                    logger.warning('Missing composites for %s %s; skipping plot', runid, cst)
                    continue

                # Check if SST is in Kelvin and convert to Celsius if necessary
                if cst_sst_mean.data.max() > 100:
                    cst_sst_mean.data -= 273.15


                cst_precip_mean = cst_precip_mean - ndjf_precip
                cst_sst_mean = cst_sst_mean - ndjf_sst
                cst_u850_mean = cst_u850_mean - ndjf_u850
                cst_v850_mean = cst_v850_mean - ndjf_v850

            skip = 2 if runid == 'obs' else 2
            u = cst_u850_mean[::skip, ::skip]
            v = cst_v850_mean[::skip, ::skip]

            # Precip wind plots
            plt.figure()
            ax = plt.axes(projection=ccrs.PlateCarree())
            figname = os.path.join(cs_index_out_dir, f'{runid}_{cst}_precip_winds850_composite.{format}')

            cmap = plt.get_cmap('BuGn') if cst == 'NDJF' else plt.get_cmap('RdBu')
            levels = np.arange(0, 32, 4) if cst == 'NDJF' else np.arange(-10, 12, 2)
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            cf = iplt.pcolormesh(cst_precip_mean, cmap=cmap, norm=norm)
            plt.colorbar(cf, orientation='vertical')

            q = iplt.quiver(u, v, scale=2, scale_units='xy')
            ax.quiverkey(q, X=132, Y=-12, U=5, label='5 m/s', labelpos='E', coordinates='data')
            plt.ylim([-10, 25])
            plt.xlim([95, 130])
            gl = ax.gridlines(draw_labels=True, alpha=0.5)
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {'size': 8}
            gl.ylabel_style = {'size': 8}
            plt.title(f'{dir_label} {cst} UV850, PRECIP')
            plt.gca().coastlines()
            save_figure(figname, bbox_inches='tight', pad_inches=0)
            plt.close()


            # SST Wind plot
            plt.figure()
            logger.info(f"Plotting SST wind composite for {runid} {cst}")
            logger.info(f"SST composite min, max value for {runid} {cst}: {cst_sst_mean.data.min()}, {cst_sst_mean.data.max()}")
            ax = plt.axes(projection=ccrs.PlateCarree())
            figname = os.path.join(cs_index_out_dir, f'{runid}_{cst}_sst_winds850_composite.{format}')
            cmap = plt.get_cmap('hot_r') if cst == 'NDJF' else plt.get_cmap('RdBu_r')
            levels = np.arange(25, 30, 0.25) if cst == 'NDJF' else np.arange(-1, 1.1, 0.1)
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            cf = iplt.pcolormesh(cst_sst_mean, cmap=cmap, norm=norm)
            plt.colorbar(cf, orientation='vertical')

            q = iplt.quiver(u, v, scale=2, scale_units='xy')
            ax.quiverkey(q, X=132, Y=-12, U=5, label='5 m/s', labelpos='E', coordinates='data')
            plt.ylim([-10, 25])
            plt.xlim([95, 130])
            gl = ax.gridlines(draw_labels=True, alpha=0.5)
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {'size': 8}
            gl.ylabel_style = {'size': 8}
            plt.title(f'{dir_label} {cst} UV850, SST')
            plt.gca().coastlines()
            save_figure(figname, bbox_inches='tight', pad_inches=0)
            plt.close()
    except Exception:
        logger.exception('Cold surge plotting failed for %s', run.get('runid', 'unknown'))
        raise


if __name__ == '__main__':
    logger.warning('This module is intended to be imported, not run directly.')
