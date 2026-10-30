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
warnings.filterwarnings(
    'ignore',
    message='.*legacy-style attribute handling.*',
    category=Warning,
)

from src import diags_eqwaves as eqw
from src import diags_precip_extremes as extreme
from src import diags_sea_cold_surges as csl1


logger = logging.getLogger(__name__)


def _time_coord(cube):
    """Return the cube's time coordinate, accepting common Iris aliases like valid_time."""
    try:
        return cube.coord('time')
    except iris.exceptions.CoordinateNotFoundError:
        for coord in cube.coords():
            if coord.name() in {'time', 'valid_time'} or coord.standard_name == 'time':
                if coord.name() != 'time':
                    coord.rename('time')
                return coord
        raise


def _has_time_coord(cube):
    """Return True if the cube exposes a usable time axis for alignment."""
    try:
        _time_coord(cube)
        return True
    except iris.exceptions.CoordinateNotFoundError:
        return False


def _time_points(cube):
    time_coord = _time_coord(cube)
    dates = list(time_coord.units.num2date(time_coord.points))
    if not dates:
        raise ValueError(f"Cube {cube.name()} has no time points")
    return dates


def check_common_time_period(cube_list):
    """Trim cubes to their shared continuous time window.

    Raises a ValueError if there are no cubes, if any remaining cube has no time
    points, or if the cubes do not overlap in time. Cubes without a time axis are
    skipped because they cannot participate in common-period alignment.
    """
    cubes = [cube for cube in cube_list if cube is not None]
    time_cubes = []
    for cube in cubes:
        if _has_time_coord(cube):
            time_cubes.append(cube)
        else:
            logger.warning(
                "Skipping cube %s in common-time alignment because it has no usable time coordinate.",
                cube.name() if cube.name() else '<unnamed>',
            )

    if not time_cubes:
        raise ValueError("No cubes with a usable time axis were provided for time alignment")

    cube_dates = [_time_points(cube) for cube in time_cubes]
    common_start_date = max(dates[0] for dates in cube_dates)
    common_end_date = min(dates[-1] for dates in cube_dates)
    logging.info(
        "Common time period for %s cubes: %s to %s",
        len(cubes),
        common_start_date.strftime('%Y-%m-%d'),
        common_end_date.strftime('%Y-%m-%d'),
    )
    if common_start_date > common_end_date:
        raise ValueError(
            "Cubes do not share a common time period: "
            f"{common_start_date} to {common_end_date}"
        )

    aligned_cubes = []
    for cube in time_cubes:
        aligned_cube = cube.extract(
            iris.Constraint(time=lambda cell: common_start_date <= cell.point <= common_end_date)
        )
        time_coord = _time_coord(aligned_cube)
        if aligned_cube is None or time_coord.points.size == 0:
            raise ValueError(f"Cube {cube.name()} could not be aligned to the common time window")
        logger.info(
            "Cube %s has %s time points after alignment.",
            aligned_cube.name(),
            len(time_coord.points),
        )
        aligned_cubes.append(aligned_cube)

    time_lengths = [len(_time_coord(cube).points) for cube in aligned_cubes]
    if not all(length == time_lengths[0] for length in time_lengths):
        raise ValueError(
            "Cubes do not have the same number of time points after alignment: "
            f"{time_lengths}"
        )

    return aligned_cubes


def _ensure_output_dirs(sea_plots_dir, mjo_plots_dir, dir_label):
    if not sea_plots_dir or not mjo_plots_dir:
        raise ValueError("sea_plots_dir and mjo_plots_dir must both be provided")

    out_plot_dir = os.path.join(sea_plots_dir, dir_label)
    out_mjo_dir = os.path.join(mjo_plots_dir, dir_label)
    os.makedirs(out_plot_dir, exist_ok=True)
    os.makedirs(out_mjo_dir, exist_ok=True)
    return out_plot_dir, out_mjo_dir


def _run_output_label(run):
    dir_label = run['dir_label'].replace(' ', '_')
    return f"{run['runid']}_{dir_label}"


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
            logger.info("Using data directory %s for %s", candidate, runid)
            return candidate

    raise FileNotFoundError(
        f"Could not locate data directory for {runid}. Tried: {', '.join(candidates)}"
    )


def _load_cube(file_path, variable_name, long_name=None, convert_units=None):
    if not os.path.exists(file_path):
        logger.warning("Missing %s file: %s", variable_name, file_path)
        return None

    try:
        cube = iris.load_cube(file_path)
    except Exception:
        logger.exception("Failed to load %s from %s", variable_name, file_path)
        return None

    try:
        _time_coord(cube)
    except iris.exceptions.CoordinateNotFoundError:
        for coord in cube.coords():
            if coord.name() in {'time', 'valid_time'} or coord.standard_name == 'time':
                if coord.name() != 'time':
                    coord.rename('time')
                break

    if long_name is not None:
        cube.long_name = long_name

    if convert_units is not None:
        try:
            cube.convert_units(convert_units)
        except Exception:
            logger.exception(
                "Loaded %s from %s but failed to convert units to %s",
                variable_name,
                file_path,
                convert_units,
            )
            return None

    return cube


def _run_file_paths(run, data_root):
    runid = run['runid']
    return {
        'precip': os.path.join(data_root, f'{runid}_PRECIP.pp.nc'),
        'olr': os.path.join(data_root, f'{runid}_OLR.pp.nc'),
        'u850': os.path.join(data_root, f'{runid}_U850.pp.nc'),
        'v850': os.path.join(data_root, f'{runid}_V850.pp.nc'),
        'sst': os.path.join(data_root, f'{runid}_SST.pp.nc'),
    }


def _log_cube_summary(cube, name):
    """Log a lightweight summary without forcing lazy Iris data to materialise."""
    try:
        core_data = cube.core_data()
        if hasattr(core_data, 'chunks') and getattr(core_data, 'chunks', None) is not None:
            logger.info(
                "Loaded %s cube with lazy data; skipping min/max materialisation for large array.",
                name.upper(),
            )
            return

        logger.info("%s min/max: %s / %s", name.upper(), core_data.min(), core_data.max())
    except Exception:
        logger.warning("Could not compute summary stats for %s cube", name, exc_info=True)


def _cube_specs(runid):
    return {
        'olr': ('toa_outgoing_longwave_flux', None),
        'u850': ('x_wind_850hPa', None),
        'v850': ('y_wind_850hPa', None),
        'precip': ('precipitation_flux', None if runid == 'obs' else 'kg m-2 day-1'),
        'sst': ('surface_temperature', None),
    }


def _load_run_cubes(run, data_root):
    runid = run['runid']
    file_paths = _run_file_paths(run, data_root)
    specs = _cube_specs(runid)

    cubes = {}
    for name, file_path in file_paths.items():
        long_name, units = specs[name]
        cube = _load_cube(file_path, name, long_name=long_name, convert_units=units)
        if cube is None:
            continue

        if not _has_time_coord(cube):
            logger.warning(
                "Skipping %s cube %s because it has no time axis and cannot be aligned to a common period.",
                name,
                file_path,
            )
            continue

        if name in {'precip', 'sst'}:
            _log_cube_summary(cube, name)
        cubes[name] = cube

    return cubes


def _align_loaded_cubes(cubes):
    ordered_names = list(cubes.keys())
    aligned = check_common_time_period(list(cubes.values()))
    return {name: cube for name, cube in zip(ordered_names, aligned)}


def _common_period_cache_dir(data_root):
    return os.path.join(data_root, 'common_time_period')


def _common_period_cache_path(cache_dir, runid, cube_name):
    return os.path.join(cache_dir, f'{runid}_{cube_name}_common_time_period.nc')


def _load_common_period_cubes(runid, data_root, cube_names):
    cache_dir = _common_period_cache_dir(data_root)
    specs = _cube_specs(runid)
    cached_cubes = {}

    for name in cube_names:
        cache_path = _common_period_cache_path(cache_dir, runid, name)
        if not os.path.exists(cache_path):
            return None

        try:
            cube = iris.load_cube(cache_path)
        except Exception:
            logger.exception("Failed to read cached common-period cube %s for %s", name, runid)
            return None

        # Cached files can be symlinks to the raw data, so re-apply the metadata
        # that downstream filenames and plots rely on.
        long_name, units = specs.get(name, (None, None))
        if long_name is not None:
            cube.long_name = long_name
        if units is not None:
            try:
                cube.convert_units(units)
            except Exception:
                logger.exception(
                    "Failed to convert cached %s cube for %s to %s", name, runid, units
                )
                return None

        cached_cubes[name] = cube

    logger.info("Loaded common-period cubes from %s", cache_dir)
    return cached_cubes


def _save_common_period_cubes(runid, data_root, cubes, source_cubes=None):
    cache_dir = _common_period_cache_dir(data_root)
    os.makedirs(cache_dir, exist_ok=True)

    for name, cube in cubes.items():
        cache_path = _common_period_cache_path(cache_dir, runid, name)
        logger.info("Writing common-period cache for %s to %s", name, cache_path)
        try:
            source_path = None
            if source_cubes is not None and name in source_cubes:
                source_path = source_cubes[name]

            if source_path and os.path.exists(source_path):
                original_cube = iris.load_cube(source_path)
                if _time_coord(cube).points.size == _time_coord(original_cube).points.size:
                    if os.path.lexists(cache_path):
                        os.unlink(cache_path)
                    os.symlink(source_path, cache_path)
                    logger.info("Linked common-period cache %s to original file %s", cache_path, source_path)
                    continue

            iris.save(cube, cache_path)
            logger.info("Saved common-period cache for %s to %s", name, cache_path)
        except Exception:
            logger.exception("Failed to write common-period cube %s for %s", name, runid)


def _prepare_common_period_cubes(runid, data_root, cubes):
    cached_cubes = _load_common_period_cubes(runid, data_root, list(cubes.keys()))
    if cached_cubes is not None:
        return cached_cubes

    source_cubes = {}
    for name, file_path in _run_file_paths({'runid': runid}, data_root).items():
        if name in cubes and os.path.exists(file_path):
            source_cubes[name] = file_path

    aligned_cubes = _align_loaded_cubes(cubes)
    _save_common_period_cubes(runid, data_root, aligned_cubes, source_cubes=source_cubes)
    return aligned_cubes


def _has_required_cubes(cubes, required_names):
    missing = [name for name in required_names if name not in cubes]
    return missing


CS_TYPES = ['NDJF', 'CS', 'CES', 'MS', 'ES']
EQW_VARNAMES = ['precipitation_flux', 'toa_outgoing_longwave_flux', 'x_wind_850hPa', 'y_wind_850hPa']
EQW_WAVENAMES = ['MJO', 'Kelvin', 'Rossby']
EQW_SEASONS = ['ndjf', 'jjas']
EXTREME_CS_TYPES = ['NDJF', 'CS', 'CES']
EXTREME_MJO_PHASES = ['NDJF', 'MJO_234', 'MJO_678']


def expected_plot_files(run, out_plot_dir, mjo_out_dir,
                        cs_level1=True, eqw_level2=True, extreme_level3=True,
                        percentile=95, format='png'):
    """Return the plot paths each enabled diagnostic block is expected to produce."""
    runid = run['runid']
    dir_label = run['dir_label'].replace(' ', '_')
    expected = []

    if cs_level1:
        expected.extend(
            os.path.join(out_plot_dir, f'{runid}_{cst}_precip_winds850_composite.{format}')
            for cst in CS_TYPES
        )
        expected.extend(
            os.path.join(out_plot_dir, f'{runid}_{cst}_sst_winds850_composite.{format}')
            for cst in CS_TYPES
        )

    if eqw_level2:
        for varname in EQW_VARNAMES:
            for season in EQW_SEASONS:
                expected.append(
                    os.path.join(out_plot_dir, f'{runid}_{varname}_ALL_{season}_stddev.{format}')
                )
            for wavename in EQW_WAVENAMES:
                for season in EQW_SEASONS:
                    expected.append(
                        os.path.join(out_plot_dir, f'{runid}_{varname}_{wavename}_{season}_stddev.{format}')
                    )
                    expected.append(
                        os.path.join(
                            out_plot_dir,
                            f'{runid}_{varname}_{wavename}_{season}_variance_percent.{format}',
                        )
                    )
                expected.append(
                    os.path.join(
                        out_plot_dir,
                        f'{dir_label}_{varname}_{wavename}_composite_hovmoller.{format}',
                    )
                )

    if extreme_level3:
        for cst in EXTREME_CS_TYPES:
            expected.append(
                os.path.join(out_plot_dir, f'{runid}_{cst}_precip_{percentile}_percentile.{format}')
            )
            if cst != 'NDJF':
                expected.append(
                    os.path.join(
                        out_plot_dir,
                        f'{runid}_{cst}_precip_{percentile}_percentile_fraction.{format}',
                    )
                )
        for mjo_phase in EXTREME_MJO_PHASES:
            expected.append(
                os.path.join(mjo_out_dir, f'{runid}_{mjo_phase}_precip_{percentile}_percentile.{format}')
            )
            if mjo_phase != 'NDJF':
                expected.append(
                    os.path.join(
                        mjo_out_dir,
                        f'{runid}_{mjo_phase}_precip_{percentile}_percentile_fraction.{format}',
                    )
                )

    return expected


def missing_plot_files(run, out_plot_dir, mjo_out_dir, **kwargs):
    """Return the expected plots that are not yet on disk for a run."""
    return [
        path
        for path in expected_plot_files(run, out_plot_dir, mjo_out_dir, **kwargs)
        if not os.path.exists(path)
    ]


def _run_cold_surge_diagnostics(run, cubes, metrics, out_plot_dir):
    required = ['u850', 'v850', 'precip', 'sst']
    missing = _has_required_cubes(cubes, required)
    if missing:
        logger.warning(
            "Skipping cold surge diagnostics for %s because these cubes are missing: %s",
            run['runid'],
            ', '.join(missing),
        )
        return

    try:
        cs_metrics = csl1.cold_surge_stats(cubes['u850'], cubes['v850'], run=run, out_plot_dir=out_plot_dir)
        if cs_metrics:
            metrics.update(cs_metrics)

        csl1.cold_surge_composites(
            [cubes['u850'], cubes['v850'], cubes['precip'], cubes['sst']],
            cstype=['NDJF', 'CS', 'CES', 'MS', 'ES'],
            run=run,
            out_plot_dir=out_plot_dir,
        )
        csl1.plot_cold_surge_composites(
            cstype=['NDJF', 'CS', 'CES', 'MS', 'ES'],
            run=run,
            out_plot_dir=out_plot_dir,
        )
    except Exception:
        logger.exception("Cold surge diagnostics failed for %s", run['runid'])


def _run_eqw_diagnostics(run, cubes, metrics, out_plot_dir):
    required = ['precip', 'olr', 'u850', 'v850']
    missing = _has_required_cubes(cubes, required)
    if missing:
        logger.warning(
            "Skipping equatorial wave diagnostics for %s because these cubes are missing: %s",
            run['runid'],
            ', '.join(missing),
        )
        return

    reference_cube = cubes['olr']
    for name in required:
        try:
            working_cube = cubes[name].copy()
            for coord_name in ['latitude', 'longitude']:
                working_cube.coord(coord_name).coord_system = reference_cube.coord(coord_name).coord_system
            working_cube = working_cube.regrid(reference_cube, iris.analysis.Linear())
            eqwave_metrics = eqw.eqwaves_compute(working_cube, out_plot_dir=out_plot_dir, run=run, compute=True)
            if eqwave_metrics:
                metrics.update(eqwave_metrics)
        except Exception:
            logger.exception("Equatorial wave diagnostics failed for %s cube %s", run['runid'], name)


def _run_extreme_precip_diagnostics(run, cubes, out_plot_dir, mjo_plots_dir):
    if 'precip' not in cubes:
        logger.warning(
            "Skipping extreme precipitation diagnostics for %s because precip is missing",
            run['runid'],
        )
        return

    try:
        extreme.precip_extremes_cs(
            cubes['precip'],
            run=run,
            percentile=95,
            compute=True,
            plot=True,
            out_plot_dir=out_plot_dir,
        )
        extreme.precip_extremes_mjo(
            cubes['precip'],
            run=run,
            percentile=95,
            compute=True,
            plot=True,
            frac_plot=True,
            mjo_index_out_dir=mjo_plots_dir,
        )
    except Exception:
        logger.exception("Extreme precipitation diagnostics failed for %s", run['runid'])


def sea_compute(varnames, control=None, expt=None, obs=None,
                cs_level1=True,
                eqw_level2=True,
                extreme_level3=True, sea_plots_dir=None, mjo_plots_dir=None):
    """
    Run the SEA diagnostics for each requested dataset.

    The function loads the available cubes, aligns them to a common time period,
    then runs each enabled diagnostic block independently. A failure in one block
    is logged and does not stop the remaining blocks from running.
    """
    logger.info("Computation starts")
    runs = [run for run in [control, expt, obs] if run is not None]

    for run in runs:
        metrics = {}

        try:
            runid = run['runid']
            output_label = _run_output_label(run)
            out_plot_dir, mjo_out_dir = _ensure_output_dirs(sea_plots_dir, mjo_plots_dir, output_label)

            missing_plots = missing_plot_files(
                run,
                out_plot_dir,
                mjo_out_dir,
                cs_level1=cs_level1,
                eqw_level2=eqw_level2,
                extreme_level3=extreme_level3,
            )
            if not missing_plots:
                logger.info(
                    "All expected plots already exist for %s in %s; skipping run without opening data",
                    runid,
                    out_plot_dir,
                )
                continue

            logger.info(
                "%s expected plots missing for %s (e.g. %s); running diagnostics",
                len(missing_plots),
                runid,
                os.path.basename(missing_plots[0]),
            )

            data_root = _resolve_data_root(run)

            logger.info("Processing run %s", runid)
            cubes = _load_run_cubes(run, data_root)
            if not cubes:
                logger.warning("No input cubes were found for %s; skipping run", runid)
                continue

            cubes = _prepare_common_period_cubes(runid, data_root, cubes)

            if cs_level1:
                _run_cold_surge_diagnostics(run, cubes, metrics, out_plot_dir)

            if eqw_level2:
                _run_eqw_diagnostics(run, cubes, metrics, out_plot_dir)

            if extreme_level3:
                _run_extreme_precip_diagnostics(run, cubes, out_plot_dir, mjo_out_dir)

            try:
                csl1.print_dict(metrics, run, out_plot_dir)
            except Exception:
                logger.exception("Failed to print metrics summary for %s", runid)

            metrics_file = os.path.join(out_plot_dir, 'metrics.csv')
            try:
                with open(metrics_file, 'w', newline='') as fh:
                    writer = csv.writer(fh)
                    for metric in metrics.items():
                        writer.writerow(metric)
                logger.info("Written metrics to %s", metrics_file)
            except Exception:
                logger.exception("Failed to write metrics CSV for %s", runid)

        except Exception:
            logger.exception("Unhandled error while processing run %s", run.get('runid', 'unknown'))
