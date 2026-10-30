import glob
import logging
import os
import sys
import uuid

import iris

logger = logging.getLogger(__name__)


def _pp2nc_regrid(combined_ppfile, ncfile):
    '''
    Regrids pp file and subsets for the tropics
    :param combined_ppfile:
    :type combined_ppfile:
    :param ncfile:
    :type ncfile:
    :return:
    :rtype:
    '''
    cube = iris.load_cube(combined_ppfile)

    #regrid to 2.5 deg
    #cube = _regrid2obs(cube)

    # tropics alone
    cube = cube.intersection(latitude=(-30, 30), longitude=(0, 360))
    print(cube)
    iris.save(cube, ncfile, netcdf_format="NETCDF3_CLASSIC")

def _pp2nc_subset(combined_ppfile, ncfile):
    '''
    Convert pp file to netcdf

    :param combined_ppfile:
    :type combined_ppfile:
    :param ncfile:
    :type ncfile:
    :return:
    :rtype:
    '''
    cube = iris.load_cube(combined_ppfile)

    # tropics alone
    cube = cube.intersection(latitude=(-30, 30), longitude=(0, 360))
    print(cube)
    iris.save(cube, ncfile, netcdf_format="NETCDF3_CLASSIC")

def _resolve_data_root(run):
    """Return the most likely data directory for a run, checking runid/label variants first."""
    base_dir = run['data_retrieve_dir']
    runid = run['runid']
    dir_label = run['dir_label'].replace(' ', '_')

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

    fallback = os.path.join(base_dir, dir_label)
    logger.info('No existing data directory found for %s; using default %s', runid, fallback)
    return fallback


def _expected_data_files(data_root, runid, varname):
    """Return the PP and NetCDF outputs expected for a variable and run."""
    combined_ppfile = os.path.join(data_root, f'{runid}_{varname}.pp')
    ncfile = combined_ppfile + '.nc'
    return combined_ppfile, ncfile


def model_data_retrieve(varnames, control=None, expt=None, netcdf=True):
    '''
    Ensure that the required model data files exist for each requested variable.

    If the NetCDF output is missing, this routine retrieves the PP data from MOOSE
    (when needed), converts it to NetCDF, and cleans up the temporary PP file.
    '''
    runs = [run for run in [control, expt] if run is not None]
    logger.info('Checking data availability for %s', [run['runid'] for run in runs])

    for run in runs:
        runid = run['runid']
        dir_label = run['dir_label']
        start_date = run['start_date']
        end_date = run['end_date']
        model = run['model']

        data_root = _resolve_data_root(run)
        os.makedirs(data_root, exist_ok=True)

        filter_tmp_file = os.path.join('/var/tmp', f'diags_filter_temp_level_{uuid.uuid4()}.txt')

        for requested_varname in varnames:
            if model == 'UM':
                filt_list = glob.glob(os.path.join('data', 'UM_Filters', requested_varname + '*.filter'))
            elif model == 'LFRic':
                filt_list = glob.glob(os.path.join('data', 'LFRic_Filters', requested_varname + '*.filter'))
            else:
                filt_list = []

            if not filt_list:
                logger.warning('No filter files found for %s in model %s', requested_varname, model)
                continue

            for filter_path in filt_list:
                ffile = os.path.basename(filter_path)
                filtered_varname = ffile.split('.')[0]
                stream = ffile.split('.')[1]
                combined_ppfile, ncfile = _expected_data_files(data_root, runid, filtered_varname)

                if os.path.exists(ncfile):
                    logger.info('%s already exists; skipping retrieval', ncfile)
                    continue

                logger.info('Preparing %s for %s (%s)', filtered_varname, runid, model)
                if os.path.exists(filter_tmp_file):
                    os.remove(filter_tmp_file)

                with open(filter_path, 'r', encoding='utf-8') as src_handle:
                    lines = src_handle.readlines()

                with open(filter_tmp_file, 'w', encoding='utf-8') as dst_handle:
                    for line in lines:
                        if line.strip() == 'begin':
                            dst_handle.write(line)
                            dst_handle.write(f'T1>={{{start_date}}}\n')
                            dst_handle.write(f'T2<={{{end_date}}}\n')
                        else:
                            dst_handle.write(line)

                if not os.path.exists(combined_ppfile):
                    logger.info('%s not found; retrieving from MASS for %s', combined_ppfile, runid)
                    moodir = f'moose:/crum/{runid}/ap{stream}.pp'
                    command = (
                        f'/usr/bin/moo select --fill-gaps-and-overwrite-smaller-files '
                        f'--condense {filter_tmp_file} {moodir} {combined_ppfile}'
                    )
                    logger.info(command)
                    status = os.system(command)
                    if status != 0:
                        raise RuntimeError(f'Retrieval failed for {combined_ppfile} (exit {status})')
                    os.system(f'./utils/sortpp {combined_ppfile}')

                if netcdf and not os.path.exists(ncfile):
                    if os.path.exists(combined_ppfile):
                        logger.info('Converting %s to NetCDF %s', combined_ppfile, ncfile)
                        _pp2nc_subset(combined_ppfile, ncfile)
                    else:
                        raise FileNotFoundError(f'Expected PP file not present after retrieval: {combined_ppfile}')

                if netcdf and os.path.exists(ncfile):
                    if os.path.exists(combined_ppfile):
                        logger.info('Removing temporary PP file %s after NetCDF conversion', combined_ppfile)
                        os.remove(combined_ppfile)
                elif not netcdf and os.path.exists(combined_ppfile):
                    logger.info('Keeping available PP file %s because netcdf=False', combined_ppfile)

