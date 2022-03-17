#!/opt/scitools/environments/default/2018_05_22-1/bin/python
import src.retrieve_SEA_data as retrieve
import do_sea as do_sea
import do_bsiso as bsiso
import mjo.do_mjo as do_mjo

if __name__ == '__main__':
    '''
    SEAPy configuration and main file.
    
    Edit the following sections and run as 
    python main.py
    
    Generate dictionary of run meta data
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

    varnames = ['U850', 'U200', 'OLR', 'V850', 'PRECIP', 'SST']

    '''
    # 1. Data download from MASS
    This is relevant for Met Office only.
    Data retrieved as daily data for the whole period for each variable
    in to separate files e.g. runid_varname.pp in to the data_retrieve_dir folder
    '''
    retrieve.model_data_retrieve(varnames, control=control, expt=expt1)

    '''
    # 2. Do MJO calculations
    MJO computation is now part of SEAPy.
     
    It has 4 sections:
    1. Mean and variance plots
    2. Wavenumber-Frequency spectra
    3. Wheeler-Hendon plots and phase composites
    4. Propagation plots and statistics
    
    Set obs=None if you do not wish to compute obs every time
    '''
    #do_mjo.mjo_compute(control=control, expt=expt1, obs=obs,
    #                  level1=True, level2=True, level3=True,
    #                  level4_prop=True)

    '''
    # 3. Do SEA computations
    This section does regional features assessment for SEAsia.
    There are 3 main sections:
    1. NDJFM mean, variance, and Cold surge types, statistics and composites
    2. Equatorial waves - Wavenumber-frequency filtering, computation of variance in each
       wave domain, ratio of wave variance to total variance, propagation composites
    3. Extreme precip statistics, % changes due to cold surge types, MJO, equatorial waves
    
    Set obs=None if you do not wish to compute obs every time
    '''
    do_sea.sea_compute(varnames, control=None, expt=None, obs=obs,
                       cs_level1=True, eqw_level2=True,
                       extreme_level3=True)

    '''
    4. BSISO computations on high res data
    This section does the computation of native data resolution (no regridding performed)
    
    Computes a simple index of NH summer ISO. 5 stages to this section
    1. filters the data to 10-90? days band
    2. compute peaks of area average time series and dates
    3. lead-lag correlations with respect to peak dates
    4. Compute extremes at different phases of ISO
    5. Generate plots of composites
    '''
    bsiso.bsiso_compute(control=control, expt=expt1, obs=None,
                        stage1_filter_variance=True,
                        stage2_iso_peaks=True,
                        stage3_iso_lag_composite=True,
                        stage4_compute_extremes=True,
                        stage5_plot_comp=True)

