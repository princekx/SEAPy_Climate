#!/opt/scitools/environments/default/2019_02_27/bin/python
import retrieve_MJO_data as rmd
import do_mjo as do_mjo

if __name__ == '__main__':
    control = {}
    control['runid'] = 'u-bz877'
    control['start_date'] = '1981/12/01'
    control['end_date'] = '2008/12/01'
    control['data_retrieve_dir'] = '/scratch/hadpx/hadgem3/data/MJOPy'

    expt = {}
    expt['runid'] = 'u-ce944'
    expt['start_date'] = '1981/12/01'
    expt['end_date'] = '2008/12/01'
    expt['data_retrieve_dir'] = '/scratch/hadpx/hadgem3/data/MJOPy'

    obs = {}
    obs['runid'] = 'obs'
    obs['start_date'] = '1989/01/01'
    obs['end_date'] = '1998/12/01'
    obs['start_year'] = 1989
    obs['end_year'] = 1998
    obs['data_retrieve_dir'] = '/project/MJO_GCSS/hadgem3/data/obs/SEAPy_data'

    print(control, expt)

    # 1. Data download from MASS
    # Check data_paths.py for data locations
    rmd.model_data_retrieve(control=control, expt=expt)

    # 2. Do MJO
    # put obs=None if you do not wish to compute obs every time
    do_mjo.mjo_compute(control=control, expt=expt, obs=None,
                       level1=False, level2=True, level3=False)
    # For Yoko, u-ai674 for control and u-bd422 for experiment.
