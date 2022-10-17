import iris
import os, sys
import numpy as np
from scipy import stats
#from tqdm import tqdm
import matplotlib.pyplot as plt
from . import mjo_utils as mu
from . import mjo_plots as mjp

class plot_contour_levels():
    def __init__(self, varname):
        # Define contour levels for plots
        if varname == 'toa_outgoing_longwave_flux':
            self.clevels = list(range(180, 280, 10))
            self.cmap = 'YlGnBu_r'
        elif varname == 'toa_outgoing_longwave_flux_filt':
            self.clevels = list(range(-20, 24, 4))
            self.cmap = 'RdBu_r'
        elif varname == 'precipitation_flux':
            self.clevels = list(range(0, 16, 1))
            self.cmap = 'YlGnBu'
        elif varname == 'precipitation_flux_filt':
            self.clevels = list(range(-5, 6, 1))
            self.cmap = 'RdBu'
        elif varname == 'x_wind_850hPa':
            self.clevels = list(range(-10, 12, 2))
            self.cmap = 'RdBu'
        elif varname == 'x_wind_850hPa_filt':
            self.clevels = list(range(-3, 4, 1))
            self.cmap = 'RdBu'
        elif varname == 'x_wind_200hPa':
            self.clevels = list(range(-20, 24, 4))
            self.cmap = 'RdBu'
        elif varname == 'x_wind_200hPa_filt':
            self.clevels = list(range(-5, 6, 1))
            self.cmap = 'RdBu'

def lagcorr(x, y, lag=None, verbose=False):
    """Compute lead-lag correlations between 2 time series.
    <x>,<y>: 1-D time series.
    <lag>: lag option, could take different forms of <lag>:
          if 0 or None, compute ordinary correlation and p-value;
          if positive integer, compute lagged correlation with lag
          upto <lag>;
          if negative integer, compute lead correlation with lead
          upto <-lag>;
          if pass in an list or tuple or array of integers, compute
          lead/lag correlations at different leads/lags.

    Note: when talking about lead/lag, uses <y> as a reference.
    Therefore positive lag means <x> lags <y> by <lag>, computation is
    done by shifting <x> to the left hand side by <lag> with respect to
    <y>.
    Similarly negative lag means <x> leads <y> by <lag>, computation is
    done by shifting <x> to the right hand side by <lag> with respect to
    <y>.

    Return <result>: a (n*2) array, with 1st column the correlation
    coefficients, 2nd column correpsonding p values.

    Currently only works for 1-D arrays.
    """
    if len(x) != len(y):
        print('Input variables in lagcorr(x,y) are of different lengths.')
        sys.exit()
    # --------Unify types of <lag>-------------
    if np.isscalar(lag):
        if abs(lag) >= len(x):
            raise Exception('Maximum lag equal or larger than array.')
        if lag < 0:
            lag = -np.arange(abs(lag) + 1)
        elif lag == 0:
            lag = [0, ]
        else:
            lag = np.arange(lag + 1)
    elif lag is None:
        lag = [0, ]
    else:
        lag = np.asarray(lag)
    # -------Loop over lags---------------------
    result = []
    if verbose:
        print('<lagcorr>: Computing lagged-correlations at lags:', lag)
    for ii in lag:
        if ii < 0:
            result.append(stats.pearsonr(y[:ii], x[-ii:]))
        elif ii == 0:
            result.append(stats.pearsonr(x, y))
        elif ii > 0:
            result.append(stats.pearsonr(y[ii:], x[:-ii]))
    result = np.asarray(result)
    return result


def find_EP_ED_ENS_STATS_OLR_3TS(var, run, out_plot_dir):
    runid = run['runid']
    loni, lati, lonf, latf = [50, -20, 180, 20]
    var = var.intersection(longitude=(loni, lonf), latitude=(lati, latf))

    var = var[51:-51]  # to avoid missing values

    lon_ens = np.linspace(-10, 10, 11)  # [-10, -7.5, -5]#, -2.5, 0, 2.5, 5, 7.5, 10]
    lat_ens = np.linspace(-5, 5, 11)  # [-5, -2.5, 0, 2.5, 5]
    npds = list(range(5, 10, 1))
    #################################
    time = var.coord('time')
    year, month, day = mu.getDates(time)
    yyyymmdd = year * 10000 + month * 100 + day

    ntime, = yyyymmdd.shape

    # sampling several times
    all_count = []
    io_count = []
    mc_count = []
    wp_count = []

    for npd in npds:
        print('%s/%s' %(npd, len(npds)))
        #for i in tqdm(range(len(lon_ens))):
        for i in range(len(lon_ens)):
            lon_en = lon_ens[i]

            for lat_en in lat_ens:
                # Sampling!!!
                # Indian Ocean
                loni, lati, lonf, latf = [70 + lon_en, -10 + lat_en, 100 + lon_en, 10 + lat_en]
                IO = var.intersection(longitude=(loni, lonf), latitude=(lati, latf))
                IO = IO.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)

                # Maritime Continent
                loni, lati, lonf, latf = [100 + lon_en, -10 + lat_en, 130 + lon_en, 10 + lat_en]
                MC = var.intersection(longitude=(loni, lonf), latitude=(lati, latf))
                MC = MC.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)

                # West Pacific
                loni, lati, lonf, latf = [130 + lon_en, -10 + lat_en, 160 + lon_en, 10 + lat_en]
                WP = var.intersection(longitude=(loni, lonf), latitude=(lati, latf))
                WP = WP.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)

                # Normalise
                IO = (IO.data - np.average(IO.data)) / np.std(IO.data)
                MC = (MC.data - np.average(MC.data)) / np.std(MC.data)
                WP = (WP.data - np.average(WP.data)) / np.std(WP.data)

                # Lead-lag correlation with IO time series
                lags = np.arange(-30, 31, 1)
                ccIO = lagcorr(IO, IO, lags)[:, 0]
                ccMC = lagcorr(IO, MC, lags)[:, 0]
                ccWP = lagcorr(IO, WP, lags)[:, 0]

                # Finding the maximum lag with respect to IO
                IO_maxlag = lags[np.where(ccIO == max(ccIO))[0]]
                MC_maxlag = lags[np.where(ccMC == max(ccMC))[0]]
                WP_maxlag = lags[np.where(ccWP == max(ccWP))[0]]

                # print IO_maxlag, MC_maxlag, WP_maxlag
                i = npd
                all_ind = []
                io_ind = []
                mc_ind = []
                wp_ind = []
                while i <= ntime - npd - 1:
                    if (month[i] >= 11 or month[i] <= 4) and \
                            (IO[i - 1] > IO[i] and IO[i + 1] > IO[i] and IO[i] <= -1.0):
                        all_ind.append(i)

                        # Check condition for MC MJO
                        if MC[i + npd] <= -1.0:
                            # Check condition for WP MJO
                            if WP[i + npd] <= -1.0:
                                wp_ind.append(i)
                                i = i + npd
                            else:
                                mc_ind.append(i)
                                i = i + npd
                        else:
                            io_ind.append(i)
                        i = i + npd
                    else:
                        i += 1

                all_count.append(len(all_ind))
                io_count.append(len(io_ind))
                mc_count.append(len(mc_ind))
                wp_count.append(len(wp_ind))

    #print(all_count)
    #print(io_count)
    #print(mc_count)
    #print(wp_count)

    print(np.percentile(all_count, 25), np.percentile(all_count, 50), np.mean(all_count), np.percentile(all_count, 75))
    print(np.percentile(io_count, 25), np.percentile(io_count, 50), np.mean(io_count), np.percentile(io_count, 75))
    print(np.percentile(mc_count, 25), np.percentile(mc_count, 50), np.mean(mc_count), np.percentile(mc_count, 75))
    print(np.percentile(wp_count, 25), np.percentile(wp_count, 50), np.mean(wp_count), np.percentile(wp_count, 75))

    file = open(os.path.join(out_plot_dir, "STATS_ENS_%s_MJO.txt" % runid), "w")
    file.write('%s %s %s %s \n' % (
    np.percentile(all_count, 25), np.percentile(all_count, 50), np.mean(all_count), np.percentile(all_count, 75)))
    file.write('%s %s %s %s \n' % (
    np.percentile(io_count, 25), np.percentile(io_count, 50), np.mean(io_count), np.percentile(io_count, 75)))
    file.write('%s %s %s %s \n' % (
    np.percentile(mc_count, 25), np.percentile(mc_count, 50), np.mean(mc_count), np.percentile(mc_count, 75)))
    file.write('%s %s %s %s \n' % (
    np.percentile(wp_count, 25), np.percentile(wp_count, 50), np.mean(wp_count), np.percentile(wp_count, 75)))
    file.close()
    print('%s written.' %file)

def find_EP_ED_dates_OLR_3TS(var, run, out_plot_dir, plotCC=False):
    runid = run['runid']
    var = var.intersection(latitude=(-30, 30), longitude=(0,360))
    var = var[51:-51] # to avoid missing values

    #################################
    time = var.coord('time')
    year, month, day = mu.getDates(time)
    yyyymmdd = year * 10000 + month * 100 + day

    ntime, = yyyymmdd.shape

    # Indian Ocean
    IO = var.intersection(latitude=(-10, 10), longitude=(70,100))
    IO = IO.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)

    # Maritime Continent
    MC = var.intersection(latitude=(-10, 10), longitude=(100,130))
    MC = MC.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)
    #comp
    # West Pacific
    WP = var.intersection(latitude=(-10, 10), longitude=(130,160))
    WP = WP.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)


    # Normalise
    IO = (IO.data - np.average(IO.data)) / np.std(IO.data)
    MC = (MC.data - np.average(MC.data)) / np.std(MC.data)
    WP = (WP.data - np.average(WP.data)) / np.std(WP.data)


    # Lead-lag correlation with IO time series
    lags = np.arange(-30, 31, 1)
    ccIO = lagcorr(IO, IO, lags)[:, 0]
    ccMC = lagcorr(IO, MC, lags)[:, 0]
    ccWP = lagcorr(IO, WP, lags)[:, 0]

    # Finding the maximum lag with respect to IO
    IO_maxlag = lags[np.where(ccIO == max(ccIO))[0]]
    MC_maxlag = lags[np.where(ccMC == max(ccMC))[0]]
    WP_maxlag = lags[np.where(ccWP == max(ccWP))[0]]

    print(IO_maxlag, MC_maxlag, WP_maxlag)

    npd = 7#MC_maxlag
    i = npd
    all_ind = []
    io_ind = []
    mc_ind = []
    wp_ind = []
    while i <= ntime - npd:
        if (month[i] >= 11 or month[i] <= 4) and \
        (IO[i - 1] > IO[i] and IO[i + 1] > IO[i] and IO[i] <= -1.0):
            all_ind.append(i)

            # Check condition for MC MJO
            if MC[i + npd] <= -1.0 :
                # Check condition for WP MJO
                if WP[i + npd] <= -1.0 :
                    wp_ind.append(i)
                    i = i + npd
                else:
                    mc_ind.append(i)
                    i = i + npd
            else:
                io_ind.append(i)
            i = i + npd
        else:
            i += 1

    all_count = len(all_ind)
    io_count = len(io_ind)
    mc_count = len(mc_ind)
    wp_count = len(wp_ind)

    print('!!!!!!!!!!Stats....!!!!!!!!!!!!!!!!!!!!')
    print('IO = %s, MC = %s, WP = %s, Total events = %s' % (io_count, mc_count, wp_count, all_count))
    print('IO Fraction = %s' % (float(io_count) * 100. / all_count))
    print('MC Fraction = %s' % (float(mc_count) * 100. / all_count))
    print('WP Fraction = %s' % (float(wp_count) * 100. / all_count))

    # Write out dates
    file = open(os.path.join(out_plot_dir, "PROP_ALL_%s_dates.txt" % runid), "w")
    for i in all_ind:
        file.write('%s\n' % yyyymmdd[i])
    file.close()
    print('%s written.' % file)

    file = open(os.path.join(out_plot_dir, "PROP_IO_%s_dates.txt" % runid), "w")
    for i in io_ind:
        file.write('%s\n' % yyyymmdd[i])
    file.close()
    print('%s written.' % file)

    file = open(os.path.join(out_plot_dir, "PROP_MC_%s_dates.txt" % runid), "w")
    for i in mc_ind:
        file.write('%s\n' % yyyymmdd[i])
    file.close()
    print('%s written.' % file)

    file = open(os.path.join(out_plot_dir, "PROP_WP_%s_dates.txt" % runid), "w")
    for i in wp_ind:
        file.write('%s\n' % yyyymmdd[i])
    file.close()
    print('%s written.' % file)

    #
    if plotCC:
        plt.plot(lags, ccIO)
        plt.plot(lags, ccMC)
        plt.plot(lags, ccWP)
        plt.legend(['Indian Ocean', 'Maritime Continent', 'Western Pacific'])
        plt.grid()
        plt.xlabel('Lags (days)')
        plt.ylabel('Correlation coefficient')
        plt.ylim([-0.7, 1.0])
        fig_name = os.path.join(out_plot_dir, 'Corr_OLR_3TS_%s.png' % runid)
        plt.savefig(fig_name)
        plt.close()
        print('Plotted %s' %fig_name)
        #plt.show()

    return io_count, mc_count, wp_count, all_count

def return_indices_of_a(a, b):
    b_set = set(b)
    return [int(i) for i, v in enumerate(a) if v in b_set]

def compute_lag_composite_3D(var, run, out_plot_dir, case='ALL'):

    # Using a efficient way of filtering.
    # Just choose the day of each case,
    # Extract a 201 day segment around the date and filter it
    # Then composite

    # Filter weights

    extract_window = 100
    lag = 30

    runid = run['runid']
    var = var.intersection(latitude=(-15, 15), longitude=(30, 200))
    var_name = var.long_name

    time = var.coord('time')
    year, month, day = mu.getDates(time)
    yyyymmdd = year * 10000 + month * 100 + day
    ntime, = yyyymmdd.shape


    filename = os.path.join(out_plot_dir, 'PROP_%s_%s_dates.txt' % (case, runid))
    dates = np.loadtxt(filename)
    ndates = len(dates)
    print('Number of dates for %s = %s' % (case, ndates))

    inds = return_indices_of_a(yyyymmdd, dates)

    if len(inds) > 3:
        n2 = 3
    comp_cube = var[inds[n2] - lag:inds[n2] + lag + 1].copy()
    comp_data = []

    #for i in tqdm(range(len(inds))):
    for i in range(len(inds)):
        ind = inds[i]
        # print('Compositing for event %s/%s' %(i+1, len(inds)))
        try:
            data_cut = var[ind - lag:ind + lag + 1].data
            if len(data_cut) == len(comp_cube.data):
                comp_data.append(var[ind - lag:ind + lag + 1].data)
            else:
                print('Date unusable for compositing. Skip!')
                pass
        except:
            print('Date unusable for compositing. Skip!')
            pass

    comp_cube.data = np.ma.mean(np.array(comp_data), axis=0)

    print(comp_cube)
    # Writing files
    comp_out_file = os.path.join(out_plot_dir, '%s_%s_comp_%s.nc' % (runid, var_name, case))
    iris.save(comp_cube, comp_out_file)
    print('%s written.' % comp_out_file)


def plot_prop_composites(var_names, run, out_plot_dir, cases=['ALL'], perc_label=None):

    runid = run['runid']
    label = run['label']
    for var_name in var_names:
        fig = plt.figure(figsize=(12, 3.75), dpi=80)
        for i, case in enumerate(cases):
            plt.subplot(1, len(cases), i + 1)
            file_name = os.path.join(out_plot_dir, '%s_%s_comp_%s.nc' % (runid, var_name, case))
            cube = iris.load_cube(file_name)
            hov = cube.intersection(latitude=(-15, 15)).collapsed('latitude', iris.analysis.MEAN)
            times = hov.coord('time').points
            times = np.arange(len(times)) - 30
            lons = hov.coord('longitude').points
            L, T = np.meshgrid(lons, times)

            levels = plot_contour_levels(var_name).clevels
            cmap = plot_contour_levels(var_name).cmap
            CS = plt.contourf(L, T, hov.data, levels=levels, cmap=cmap, extend='both')
            plt.title('%s %s events \n%s' % (label, case, var_name))
            plt.xlabel('Longitude (degrees east)')
            if i == 0:
                plt.ylabel('Lag/lead (days)')
            plt.grid()
            if case == cases[-1]:
                pos = plt.gca().get_position()
                pos2 = fig.add_axes([pos.x0 + pos.width + 0.02, pos.y0, pos.width / 10.0, pos.height])
                plt.colorbar(CS, cax=pos2)

        plt_out_file = os.path.join(out_plot_dir, 'Plot_%s_%s_comp.png' % (runid, var_name))
        plt.savefig(plt_out_file)
        plt.close()
        print('Plotted %s .' % plt_out_file)

def compute_composites_forSFCData(vars, run, out_plot_dir):
    runid = run['runid']
    label = run['label']
    for var in vars:
        # ALL
        all_comp = compute_lag_composite_3D(var, run, out_plot_dir, case='ALL')

        # IO
        io_comp = compute_lag_composite_3D(var, run, out_plot_dir, case='IO')

        # MC
        mc_comp = compute_lag_composite_3D(var, run, out_plot_dir, case='MC')

        # WP
        wp_comp = compute_lag_composite_3D(var, run, out_plot_dir, case='WP')


def diagnos_level4_prop(outgoing_longwave_cubes, u_wind_850_cubes, u_wind_200_cubes,
                precip_cubes, run, out_plot_dir=None):
    runid = run['runid']
    label = run['label']
    varname = outgoing_longwave_cubes.long_name

    # filter OLR data
    filt_olr_cube = mu.lanczos_filter(outgoing_longwave_cubes, window=101, low=1. / 100., high=1. / 20.)

    # compute the dates of the propagating events
    nIO, nMC, nWP, nALL = find_EP_ED_dates_OLR_3TS(filt_olr_cube, run, out_plot_dir)

    # composites for the whole data unfiltered
    compute_composites_forSFCData([outgoing_longwave_cubes, u_wind_850_cubes, u_wind_200_cubes, precip_cubes],
                                  run, out_plot_dir)

    var_names = [var.long_name for var in [outgoing_longwave_cubes, u_wind_850_cubes, u_wind_200_cubes, precip_cubes]]
    plot_prop_composites(var_names, run, out_plot_dir, cases=['ALL', 'IO', 'MC', 'WP'])

    # composites for filtered cubes
    for cube in [outgoing_longwave_cubes, u_wind_850_cubes, u_wind_200_cubes, precip_cubes]:
        cube_filt = mu.lanczos_filter(cube, window=101, low=1. / 100., high=1. / 20.)
        cube_filt.long_name = cube.long_name+'_filt'

        compute_composites_forSFCData([cube_filt], run, out_plot_dir)
        # Plotting filtered data composites
        plot_prop_composites([cube_filt.long_name], run, out_plot_dir, cases=['ALL', 'IO', 'MC', 'WP'])

    return
