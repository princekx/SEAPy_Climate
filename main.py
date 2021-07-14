#!/opt/scitools/environments/default/2018_05_22-1/bin/python
import sys
import src.retrieve_SEA_data as retrieve
import do_sea as do_sea
import mjo.retrieve_MJO_data as rmd
import mjo.do_mjo as do_mjo

if __name__ == '__main__':
    control = {}
    control['runid'] = 'u-bw062'
    control['start_date'] = '1981/12/01'
    control['end_date'] = '2008/12/01'
    control['data_retrieve_dir'] = '/scratch/hadpx/hadgem3/data/SEAPy'

    #control['start_date'] = '2022/01/01'
    #control['end_date'] = '2041/12/01'
    #control['data_retrieve_dir'] = '/scratch/hadpx/hadgem3/data/SEAPy'

    expt = {}
    expt['runid'] = 'u-by782'
    expt['start_date'] = '1981/12/01'
    expt['end_date'] = '2008/12/01'
    expt['data_retrieve_dir'] = '/scratch/hadpx/hadgem3/data/SEAPy'

    obs = {}
    obs['runid'] = 'obs'
    obs['start_date'] = '1989/01/01'
    obs['end_date'] = '1998/12/01'
    obs['data_retrieve_dir'] = '/project/MJO_GCSS/hadgem3/data/obs/SEAPy_data'

    print(control, expt)
    varnames = ['U850', 'U200', 'OLR', 'V850', 'PRECIP', 'SST']

    # 1. Data download from MASS
    #retrieve.model_data_retrieve(varnames, control=control, expt=expt)

    # 2. Do MJO
    # put obs=None if you do not wish to compute obs every time
    #do_mjo.mjo_compute(control=control, expt=expt, obs=None,
    #                   level1=False, level2=False, level3=False)

    # 3. Do SEA computations
    # put obs=None if you do not wish to compute obs every time
    #do_sea.sea_compute(varnames, control=control, expt=expt, obs=obs,
    #                   cs_level1=True, eqw_level2=True, level3=False)
    do_sea.sea_compute(varnames, control=None, expt=None, obs=obs,
                       cs_level1=False, eqw_level2=False, mjo_level3=False,
                       extreme_level4=True)
