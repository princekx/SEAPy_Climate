import os
import iris
import iris
import numpy as np
from iris.coord_categorisation import add_year, add_month_number, add_day_of_month, add_day_of_year
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import iris.analysis.calculus as calculus
#from tqdm import tqdm
import datetime
from src import data_paths


def derivative(cube, axisname):
    R = 6378388.  # Radius of the earth
    deg2rad = 0.0174533  # pi/180.
    dcube = cube.copy()

    coord_names = np.array([c.var_name for c in cube.coords()])
    # print(coord_names)

    if axisname == 'latitude':
        lats = cube.coord('latitude').points
        axis_index = np.where(coord_names == 'latitude')[0][0]
        dlat = np.diff(lats) * deg2rad  # convert to radians
        dy = R * np.sin(dlat)  # constant at this latitude
        dcube = calculus.differentiate(cube, 'latitude')
        dcube /= iris.util.broadcast_to_shape(dy, dcube.shape, (axis_index,))

    if axisname == 'longitude':
        lats = cube.coord('latitude').points
        lons = cube.coord('longitude').points
        axis_index = np.where(coord_names == 'latitude')[0][0]
        dlon = (lons[1] - lons[0]) * deg2rad  # convert to radians
        dx = np.array([R * np.cos(deg2rad * lat) * dlon for lat in lats])
        dcube = calculus.differentiate(cube, 'longitude')
        dcube /= iris.util.broadcast_to_shape(dx, dcube.shape, (axis_index,))
    return dcube


def compute_vort_div(ucube, vcube, vort_file=None, div_file=None):
    print('computing vorticity and divergence...')
    if not os.path.exists(div_file):
        div = derivative(ucube, 'longitude').regrid(ucube, iris.analysis.Linear())
        div += derivative(vcube, 'latitude').regrid(ucube, iris.analysis.Linear())
        div.var_name = 'divergence'
        iris.save(div, div_file)
        print('Written %s' % div_file)
    else:
        print('File found. Skipping:  %s' % div_file)

    if not os.path.exists(vort_file):
        vort = derivative(vcube, 'longitude').regrid(ucube, iris.analysis.Linear())
        vort -= derivative(ucube, 'latitude').regrid(ucube, iris.analysis.Linear())
        vort.var_name = 'vorticity'
        iris.save(vort, vort_file)
        print('Written %s' % vort_file)
    else:
        print('File found. Skipping:  %s' % vort_file)


def prepare_calendar(cube):
    # Setting up the dates on data
    if not cube.coords('year'):
        iris.coord_categorisation.add_year(cube, 'time', name='year')
    if not cube.coords('month_number'):
        iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')
    if not cube.coords('day_of_month'):
        iris.coord_categorisation.add_day_of_month(cube, 'time', name='day_of_month')
    if not cube.coords('day_of_year'):
        iris.coord_categorisation.add_day_of_year(cube, 'time', name='day_of_year')
    return cube


def create_dates_df(cube):
    cube = prepare_calendar(cube)
    cube_dates_df = pd.concat([pd.DataFrame(cube.coord('year').points),
                               pd.DataFrame(cube.coord('month_number').points),
                               pd.DataFrame(cube.coord('day_of_month').points)], axis=1)
    cube_dates_df.columns = ['year', 'month', 'day']
    return cube_dates_df


def mean_var_season(cube, varname=None, runid=None, varlabel=None,
                    season_months=[5, 6, 7, 8, 9], writeout=True):
    '''
    Mean and variance for a given season

    :param cube:
    :type cube:
    :param varname:
    :type varname:
    :param runid:
    :type runid:
    :param varlabel:
    :type varlabel:
    :param season_months:
    :type season_months:
    :param writeout:
    :type writeout:
    :return:
    :rtype:
    '''
    cube_dates_df = create_dates_df(cube)
    season_df = cube_dates_df.loc[(cube_dates_df['month'].isin(season_months))]

    indices = [i for i in season_df.index]
    seasonal_variance = cube[indices].collapsed('time', iris.analysis.VARIANCE)
    seasonal_mean = cube[indices].collapsed('time', iris.analysis.MEAN)

    if writeout:
        cs_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
        writeout_file = os.path.join(cs_index_out_dir,
                                     '%s_MJJAS_%s_mean.nc' % (runid, varname))
        iris.save(seasonal_mean, writeout_file)
        print('Written %s' % writeout_file)

        writeout_file = os.path.join(cs_index_out_dir,
                                     '%s_MJJAS_%s_variance.nc' % (runid, varname))
        iris.save(seasonal_variance, writeout_file)
        print('Written %s' % writeout_file)


def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    """
    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n - 1:0:-1] = firstfactor * sigma
    w[n + 1:-1] = firstfactor * sigma
    return w[1:-1]


def lanczos_filter(cube, window=101, low=1. / 100., high=1. / 20.):
    '''
    Driver for a band-pass filter

    :param cube:
    :type cube:
    :param window:
    :type window:
    :param low:
    :type low:
    :param high:
    :type high:
    :return:
    :rtype:
    '''
    wgts_high = low_pass_weights(window, high)
    wgts_low = low_pass_weights(window, low)

    # remove long term mean
    cube = cube - cube.collapsed('time', iris.analysis.MEAN)

    cube_filt = cube.copy()
    cube_filt.data = np.apply_along_axis(lambda m: (np.convolve(m, wgts_high, mode='same') -
                                                    np.convolve(m, wgts_low, mode='same')), axis=0, arr=cube.data)

    # for masked data
    cube_filt.data = np.ma.masked_array(cube_filt.data, cube.data.mask)

    return cube_filt


def return_indices_of_a(a, b):
    '''
    Returns the indices of A that matches elements in B

    :param a:
    :type a:
    :param b:
    :type b:
    :return:
    :rtype:
    '''
    b_set = set(b)
    return [int(i) for i, v in enumerate(a) if v in b_set]


def return_indices_of_adf(adf, bdf):
    '''
    Returns the indices of DataFrame A that matches elements in DataFrame B

    :param adf:
    :type adf:
    :param bdf:
    :type bdf:
    :return:
    :rtype:
    '''
    a = [year * 10000 + month * 100 + day for year, month, day in zip(adf.year, adf.month, adf.day)]
    b = [year * 10000 + month * 100 + day for year, month, day in zip(bdf.year, bdf.month, bdf.day)]
    b_set = set(b)
    return [int(i) for i, v in enumerate(a) if v in b_set]


def make_dates_from_df(df):
    return [year * 10000 + month * 100 + day for year, month, day in zip(df.year, df.month, df.day)]


def make_composite(cube, peak_dates_df, runid=None, varname=None, lag=20, write_out=True):
    '''
    Composite a cube with lead-lag around a set of peaks
    :param cube:
    :type cube: Iris cube
    :param peak_dates_df:
    :type peak_dates_df: Pandas DataFrame
    :param lag: compute from -lag:+lag
    :type lag: days, integer
    :return:
    :rtype: iris cube
    '''
    peaks_dates = make_dates_from_df(peak_dates_df)
    cube_dates = make_dates_from_df(create_dates_df(cube))

    inds = return_indices_of_a(cube_dates, peaks_dates)
    # Get a mid value of the number of indices
    if len(inds) > 3:
        n2 = 5
    comp_cube = cube[inds[n2] - lag:inds[n2] + lag + 1].copy()
    comp_data = []
    #for i in tqdm(range(len(inds))):
    for i in range(len(inds)):
        ind = inds[i]
        print('Compositing %s for event %s/%s' %(varname, i+1, len(inds)))
        try:
            data_cut = cube[ind - lag:ind + lag + 1].data
            if len(data_cut) == len(comp_cube.data):
                comp_data.append(cube[ind - lag:ind + lag + 1].data)
            else:
                print('Date unusable for compositing. Skip!')
                pass
        except:
            print('Date unusable for compositing. Skip!')
            pass

    comp_cube.data = np.ma.mean(np.array(comp_data), axis=0)

    if write_out:
        out_plot_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
        print(out_plot_dir)
        writeout_file = os.path.join(out_plot_dir,
                                     '%s_MJJAS_ISO_lead_lag_%s_composite.nc' % (runid, varname))
        iris.save(comp_cube, writeout_file)
        print('Written %s' % writeout_file)


def indices_around_peaks(inds, days_around_peak=3):
    '''
    Indices of time series and 3 days around the peak

    :param inds:
    :type inds:
    :param days_around_peak:
    :type days_around_peak:
    :return:
    :rtype:
    '''
    inds_around_peak = []
    for ind in inds:
        inds_around_peak.extend(list(np.arange(ind - days_around_peak, ind + days_around_peak)))
    # remove any negative values
    inds_around_peak = [item for item in inds_around_peak if item >= 0]
    return inds_around_peak


def compute_percentiles(cube, indices, percent=95):
    '''
    Percentile values for given time indices

    :param cube:
    :type cube:
    :param indices:
    :type indices:
    :param percent:
    :type percent:
    :return:
    :rtype:
    '''
    # compute percentiles
    percentiles = cube[indices].collapsed('time', iris.analysis.PERCENTILE, percent=percent)
    return percentiles


def make_hovmoller(cube, lon_range=(100, 120), lat_range=(-10, 25), average_along='longitude'):
    '''
    Collapse cube for given latitude-longitude range and along a given axis coordinate

    :param cube:
    :type cube:
    :param lon_range:
    :type lon_range:
    :param lat_range:
    :type lat_range:
    :param average_along:
    :type average_along:
    :return:
    :rtype:
    '''
    return cube.intersection(longitude=lon_range, latitude=lat_range).collapsed(average_along, iris.analysis.MEAN)


def plot_composite(runid, lat_range=(-10, 25), lon_range=(100, 120)):
    '''
    Plot composites

    :param runid:
    :type runid:
    :param lat_range:
    :type lat_range:
    :param lon_range:
    :type lon_range:
    :return:
    :rtype:
    '''

    out_plot_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
    print(out_plot_dir)
    precip_comp_file = os.path.join(out_plot_dir,
                                    '%s_MJJAS_ISO_lead_lag_%s_composite.nc' % (runid, 'precipitation_flux'))
    olr_comp_file = os.path.join(out_plot_dir,
                                    '%s_MJJAS_ISO_lead_lag_%s_composite.nc' % (runid, 'toa_outgoing_longwave_flux'))
    u850_comp_file = os.path.join(out_plot_dir,
                                  '%s_MJJAS_ISO_lead_lag_%s_composite.nc' % (runid, 'x_wind_850'))
    v850_comp_file = os.path.join(out_plot_dir,
                                  '%s_MJJAS_ISO_lead_lag_%s_composite.nc' % (runid, 'y_wind_850'))
    precip_comp_filt_file = os.path.join(out_plot_dir,
                                         '%s_MJJAS_ISO_lead_lag_%s_composite.nc' % (runid, 'precipitation_flux_filt'))
    olr_comp_filt_file = os.path.join(out_plot_dir,
                                 '%s_MJJAS_ISO_lead_lag_%s_composite.nc' % (runid, 'toa_outgoing_longwave_flux_filt'))
    sst_comp_filt_file = os.path.join(out_plot_dir,
                                      '%s_MJJAS_ISO_lead_lag_%s_composite.nc' % (runid, 'sst_filt'))
    vort850_comp_filt_file = os.path.join(out_plot_dir,
                                          '%s_MJJAS_ISO_lead_lag_%s_composite.nc' % (runid, 'vorticity_850_filt'))
    div850_comp_filt_file = os.path.join(out_plot_dir,
                                         '%s_MJJAS_ISO_lead_lag_%s_composite.nc' % (runid, 'divergence_850_filt'))

    print(precip_comp_file, olr_comp_file, u850_comp_file, v850_comp_file)
    precip_comp = iris.load_cube(precip_comp_file)
    olr_comp = iris.load_cube(olr_comp_file)
    u850_comp = iris.load_cube(u850_comp_file)
    v850_comp = iris.load_cube(v850_comp_file)

    # U and v on the same grid
    v850_comp = v850_comp.regrid(u850_comp, iris.analysis.Linear())

    sst_filt_comp = iris.load_cube(sst_comp_filt_file)
    precip_filt_comp = iris.load_cube(precip_comp_filt_file)
    olr_filt_comp = iris.load_cube(olr_comp_filt_file)
    vort850_filt_comp = iris.load_cube(vort850_comp_filt_file)
    div850_filt_comp = iris.load_cube(div850_comp_filt_file)

    precip_comp_hov = make_hovmoller(precip_comp, lon_range=lon_range, lat_range=lat_range,
                                     average_along='longitude')
    olr_comp_hov = make_hovmoller(olr_comp, lon_range=lon_range, lat_range=lat_range,
                                     average_along='longitude')
    u850_comp_hov = make_hovmoller(u850_comp, lon_range=lon_range, lat_range=lat_range,
                                   average_along='longitude')
    v850_comp_hov = make_hovmoller(v850_comp, lon_range=lon_range, lat_range=lat_range,
                                   average_along='longitude')
    precip_comp_filt_hov = make_hovmoller(precip_filt_comp, lon_range=lon_range, lat_range=lat_range,
                                          average_along='longitude')
    olr_comp_filt_hov = make_hovmoller(olr_filt_comp, lon_range=lon_range, lat_range=lat_range,
                                          average_along='longitude')
    sst_comp_filt_hov = make_hovmoller(sst_filt_comp, lon_range=lon_range, lat_range=lat_range,
                                       average_along='longitude')
    vort850_comp_filt_hov = make_hovmoller(vort850_filt_comp, lon_range=lon_range, lat_range=lat_range,
                                           average_along='longitude')
    div850_comp_filt_hov = make_hovmoller(div850_filt_comp, lon_range=lon_range, lat_range=lat_range,
                                          average_along='longitude')
    ntime, nlat = precip_comp_hov.shape
    lag = int(ntime / 2)

    fig_name = os.path.join(out_plot_dir, '%s_MJJAS_precip_winds_sst_leadlag_comp.png' % runid)
    fig = plt.figure(1, figsize=(12, 7), dpi=100)
    plt.subplot(121)
    times = np.arange(-lag, lag + 1, 1)
    lats = precip_comp_hov.coord('latitude').points
    T, L = np.meshgrid(times, lats)
    CS = plt.contourf(T, L, precip_comp_hov.data.T, levels=np.arange(2, 42, 2), cmap='YlGnBu', extend='max')
    plt.colorbar(CS)
    # qplt.contour(u850_comp_filt_hov, levels=np.linspace(-5, 5, 11), colors='k')
    t = times
    y = v850_comp_hov.coord('latitude').points

    if (y[0] > y[1]):
        # need to reverse latitude dimension
        print('Reversing latitude')
        u850_comp_hov = iris.util.reverse(u850_comp_hov, 'latitude')
        v850_comp_hov = iris.util.reverse(v850_comp_hov, 'latitude')
        y = v850_comp_hov.coord('latitude').points

    # y values need to be of equally spaced
    dy = y[1]-y[0]
    y = np.array([y[0]+i*dy for i in range(len(y))])
    #print(y)

    u = u850_comp_hov.data.T
    v = v850_comp_hov.data.T
    speed = np.sqrt(u ** 2 + v ** 2)
    lw = speed / 5  # an adhoc scaling for plotting

    # Bounds for streamplot
    bounds = np.arange(0, 4, 1)
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    cl = plt.gca().streamplot(t, y, u, v, density=(2., 1.5), maxlength=1, cmap='GnBu', norm=norm, color=speed,
                              linewidth=lw)

    #plt.gca().quiver(t, y, u, v, cmap='GnBu', alpha=0.5)
    plt.ylim([-10, 25])
    plt.ylabel('Latitude')
    plt.xlabel('Lead/Lag (days)')
    plt.title('%s Precip., 850hPa winds [95-110E]' % runid)

    plt.subplot(122)
    lats = precip_comp_filt_hov.coord('latitude').points
    T, L = np.meshgrid(times, lats)
    CS = plt.contourf(T, L, precip_comp_filt_hov.data.T, levels=np.linspace(-5, 5, 11), cmap='RdBu', extend='both')
    plt.colorbar(CS)

    lats = sst_comp_filt_hov.coord('latitude').points
    T, L = np.meshgrid(times, lats)
    csc = plt.contour(T, L, sst_comp_filt_hov.data.T, levels=np.arange(-0.3, 0.4, 0.05), colors='k')
    plt.gca().clabel(csc, inline=True, fontsize=10, fmt='%1.1f')
    plt.contour(T, L, sst_comp_filt_hov.data.T, levels=[0], colors='k', linewidths=[3])
    plt.ylim([-10, 25])
    plt.ylabel('Latitude')
    plt.xlabel('Lead/Lag (days)')
    plt.title('%s Filtered precip., SST [95-110E]' % runid)
    plt.savefig(fig_name)
    print('%s plotted.' % fig_name)
    plt.close()

    # Second plot
    fig_name = os.path.join(out_plot_dir, '%s_MJJAS_precip_winds_vort_div850_leadlag_comp.png' % runid)
    fig = plt.figure(2, figsize=(12, 7), dpi=100)
    plt.subplot(121)
    times = np.arange(-lag, lag + 1, 1)
    lats = precip_comp_filt_hov.coord('latitude').points
    T, L = np.meshgrid(times, lats)
    CS = plt.contourf(T, L, precip_comp_filt_hov.data.T, levels=np.linspace(-5, 5, 11), cmap='RdBu', extend='both')
    plt.colorbar(CS)

    lats = vort850_comp_filt_hov.coord('latitude').points
    T, L = np.meshgrid(times, lats)
    CS1 = plt.contour(T, L, vort850_comp_filt_hov.data.T * 1e5, levels=np.linspace(-5, 5, 11), colors='k')
    plt.gca().clabel(CS1, inline=True, fontsize=10, fmt='%1.1f')
    plt.contour(T, L, vort850_comp_filt_hov.data.T * 1e5, levels=[0], colors='k', linewidths=[3])

    plt.ylabel('Latitude')
    plt.xlabel('Lead/Lag (days)')
    plt.title('Precip., 850hPa vorticity (*1e5) [95-110E]')

    plt.subplot(122)
    times = np.arange(-lag, lag + 1, 1)
    lats = precip_comp_filt_hov.coord('latitude').points
    T, L = np.meshgrid(times, lats)
    CS = plt.contourf(T, L, precip_comp_filt_hov.data.T, levels=np.linspace(-10, 10, 11), cmap='RdBu', extend='both')
    plt.colorbar(CS)

    lats = div850_comp_filt_hov.coord('latitude').points
    T, L = np.meshgrid(times, lats)
    CS1 = plt.contour(T, L, div850_comp_filt_hov.data.T * 1e5, levels=np.linspace(-1, 1, 11), colors='k')
    plt.gca().clabel(CS1, inline=True, fontsize=10, fmt='%1.1f')
    plt.contour(T, L, div850_comp_filt_hov.data.T * 1e5, levels=[0], colors='k', linewidths=[3])

    plt.ylabel('Latitude')
    plt.xlabel('Lead/Lag (days)')
    plt.title('Precip., 850hPa divergence (*1e5) [95-110E]')
    plt.savefig(fig_name)
    plt.close()
    print('%s plotted.' % fig_name)
