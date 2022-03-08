#!/opt/scitools/environments/default/2018_05_22-1/bin/python
import sys
import src.retrieve_SEA_data as retrieve
import do_sea as do_sea
import do_bsiso as bsiso
import mjo.do_mjo as do_mjo

if __name__ == '__main__':
    '''
    SEAPy configuration
    Generate 
    
    '''
    control = {}
    control['runid'] = 'u-bu357'
    control['label'] = 'GA8GL9'
    control['start_date'] = '1981/12/01'
    control['end_date'] = '2008/12/01'
    control['data_retrieve_dir'] = '/scratch/hadpx/hadgem3/data/SEAPy'

    expt = {}
    expt['runid'] = 'u-ci336'
    expt['label'] = 'GA8GL9_577.6'
    expt['start_date'] = '1981/12/01'
    expt['end_date'] = '2008/12/01'
    expt['data_retrieve_dir'] = '/scratch/hadpx/hadgem3/data/SEAPy'

    expt1 = {}
    expt1['runid'] = 'u-ch221'
    expt1['label'] = 'CoMv8p1'
    expt1['start_date'] = '1981/12/01'
    expt1['end_date'] = '2008/12/01'
    expt1['data_retrieve_dir'] = '/scratch/hadpx/hadgem3/data/SEAPy'

    obs = {}
    obs['runid'] = 'obs'
    obs['label'] = 'ERAInt/TRMM'
    obs['start_date'] = '1989/01/01'
    obs['end_date'] = '1998/12/01'
    obs['data_retrieve_dir'] = '/project/MJO_GCSS/hadgem3/data/obs/SEAPy_data'

    print(control, expt)
    varnames = ['U850', 'U200', 'OLR', 'V850', 'PRECIP', 'SST']

    # 1. Data download from MASS
    #retrieve.model_data_retrieve(varnames, control=control, expt=expt)

    # 2. Do MJO calculations
    # put obs=None if you do not wish to compute obs every time
    #do_mjo.mjo_compute(control=control, expt=expt1, obs=obs,
    #                  level1=False, level2=False, level3=True,
    #                  level4_prop=True)

    # 3. Do SEA computations
    # put obs=None if you do not wish to compute obs every time
    do_sea.sea_compute(varnames, control=None, expt=expt, obs=None,
                       cs_level1=False, eqw_level2=False,
                       extreme_level4=True)

    # BSISO computations on high res data
    #bsiso.bsiso_compute(control=None, expt=expt, obs=None,
    #                    stage1_filter_variance=False,
    #                    stage2_iso_peaks=False,
    #                    stage3_iso_lag_composite=True,
    #                    stage4_compute_extremes=False,
    #                    stage5_plot_comp=True)

