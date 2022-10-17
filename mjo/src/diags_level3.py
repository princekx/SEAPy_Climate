import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import iris
import iris.plot as iplt
import os, sys
import pandas as pd
import datetime

from .mycolormaps import getcolors
from . import mjo_utils as mu
from . import mjo_plots as mjp

UNDEF = -999.


def find_array_match_inds(long_array, short_array):
    return long_array.searchsorted(short_array)


def return_indices_of_a(a, b):
    '''
    Function to return indices of elements in A that matches values in B
    :param a: list a
    :type a:
    :param b: list b
    :type b:
    :return: indices of elements in a that matches values in b
    :rtype:
    '''
    b_set = set(b)
    return [i for i, v in enumerate(a) if v in b_set]


def _makecube_lcorr(var, lags, lons):
    var_cube = iris.cube.Cube(var)
    var_cube.rename('Lag_correlation')
    lags_coord = iris.coords.DimCoord(lags, long_name='lead_lag')
    lags_coord.guess_bounds()
    lons_coord = iris.coords.DimCoord(lons, long_name='longitude')
    lons_coord.guess_bounds()
    var_cube.add_dim_coord(lags_coord, 0)
    var_cube.add_dim_coord(lons_coord, 1)
    return var_cube


def diagnos_level3(olr_cube, x_wind_850_cube, x_wind_200_cube, precip_cube,
                   runid, label, out_dir):
    print('Starting Level 3 diagnostics...')

    # RMM - Real time multivariate MJO Index
    rmmfile = os.path.join(out_dir, 'RMMs_' + runid + '.txt')

    # Time-filter cube data
    filtered_olr_cube = mu.Filter(olr_cube)
    filtered_x_wind_850_cube = mu.Filter(x_wind_850_cube)
    filtered_x_wind_200_cube = mu.Filter(x_wind_200_cube)
    filtered_precip_cube = mu.Filter(precip_cube)

    print(filtered_precip_cube)
    ### Lead-Lag Correlation Plot
    for cube in [filtered_olr_cube,
                 filtered_x_wind_850_cube,
                 filtered_x_wind_200_cube,
                 filtered_precip_cube]:
        varname = cube.long_name

        # if varname == 'x_wind':
        #    assert len(cube.coord('pressure').points) == 1
        #    pressure_level = cube.coord('pressure').points[0]
        #    if pressure_level == 850:
        #        varname = 'x_wind_850hPa'
        #    if pressure_level == 200:
        #        varname = 'x_wind_200hPa'

        # Extract area-average timeseries
        lon1, lat1, lon2, lat2 = [80, -10, 100, 10]
        reference_time_series = mu.AreaAverage(cube, [lon1, lat1, lon2, lat2])

        # lat average
        lon1, lat1, lon2, lat2 = [40, -10, 180, 10]
        cube = cube.extract(mu.region([lon1, lat1, lon2, lat2]))
        cube = cube.collapsed(['latitude'], iris.analysis.MEAN)
        lead_lag_corr, lags, longitudes = mu.LeadLagCorr(reference_time_series, cube)

        title = ' '.join(
            [label, "OLR [80-100E, 10S-10N] vs", varname, "Lead-Lag Correlation"]
        )

        out_name = runid + "_" + varname + "_LeadLagCorr"
        figname = os.path.join(out_dir, "%s.png" % out_name)
        mjp.HovTimeLon(lead_lag_corr, lags, longitudes, levels=np.arange(-1, 1.1, 0.1),
                       title=title, figname=figname)

        lead_lag_corr_cube = _makecube_lcorr(lead_lag_corr, lags, longitudes)
        ncname = os.path.join(out_dir, "%s.nc" % out_name)
        iris.save(lead_lag_corr_cube, ncname)

    ### Wheeler Hendon Plot
    lon1, lat1, lon2, lat2 = [0, -15, 360, 15]
    out_name = runid + "_WH04"
    figname = os.path.join(out_dir, "%s.png" % out_name)
    nwgt = 101
    n_2 = int((nwgt - 1) / 2)

    # Create large cube to hold all three variables
    cube = filtered_olr_cube

    # lat average
    cube = cube.extract(mu.region([lon1, lat1, lon2, lat2]))
    cube = cube.collapsed(['latitude'], iris.analysis.MEAN)

    time_coord = cube.coord('time')
    year, month, day = mu.getDates(time_coord)

    # Large array to hold all the variables
    ntime, mlon = cube.shape
    cdata = np.zeros((ntime, 3 * mlon))

    # Merdional mean
    for i, cube in enumerate([filtered_olr_cube,
                              filtered_x_wind_850_cube,
                              filtered_x_wind_200_cube]):
        varname = cube.long_name

        # if varname == 'x_wind':
        #    assert len(cube.coord('pressure').points) == 1
        #    pressure_level = cube.coord('pressure').points[0]
        #    if pressure_level == 850:
        #        varname = 'x_wind_850hPa'
        #    if pressure_level == 200:
        #       varname = 'x_wind_200hPa'

        # lat average
        cube = cube.extract(mu.region([lon1, lat1, lon2, lat2]))
        cube = cube.collapsed(['latitude'], iris.analysis.MEAN)

        print(cube)

        # Compute the temporal variance
        variance = np.var(cube.data, axis=0)

        # Compute the zonal mean of the temporal variance
        zonal_mean_variance = np.average(variance)

        # Normalize by standard deviation
        cube = cube / np.sqrt(zonal_mean_variance)  # time, lon

        # Combine the normalized data into one variable
        print(i * mlon, (i + 1) * mlon)
        cdata[:, i * mlon:(i + 1) * mlon] = cube.data

    # Compute Combined EOF
    # Memory failure. So reading Wheeler & Hendon 2004 EOFs directly
    # covariance matrix
    # var   = np.dot(cdata,cdata.T)
    # d , v = np.linalg.eig(var)
    # print d.shape,  v.shape

    # Read Observed EOF
    EOFs_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'WH04_EOFs.dat')
    eofs = np.loadtxt(EOFs_file)

    # DEBUG
    # import pickle
    # fileObject = open('Pickled.pkl', 'wb')
    # pickle.dump([cdata, eofs], fileObject)
    # fileObject.close()

    # Compute RMM1 and RMM2 (the first two normalized PCs)
    pcs = np.dot(cdata, eofs)

    print(pcs)
    ntime, neofs = pcs.shape

    # Now normalize (by EOF-calculated s.d.) the newly calculated PCs
    for n in np.arange(neofs):
        pcs[:, n] = (pcs[:, n] - np.average(pcs[:, n])) / np.std(pcs[:, n])

    inds = [i for i in range(ntime) if cdata[i, 10] == UNDEF]

    # Compute amplitude
    amp = np.sqrt(pcs[:, 0] ** 2 + pcs[:, 1] ** 2)

    # Compute phase
    pha = mu.RMM_phases(pcs)

    pcs[:n_2 + 1, :] = UNDEF
    pcs[ntime - n_2:, :] = UNDEF
    amp[:n_2 + 1] = UNDEF
    amp[ntime - n_2:] = UNDEF
    pha[:n_2 + 1] = UNDEF
    pha[ntime - n_2:] = UNDEF

    # open file to write RMMs into
    with open(rmmfile, 'w') as f:
        for it in np.arange(ntime):
            f.write(' '.join(
                [str(year[it]),
                 str(month[it]),
                 str(day[it]),
                 str(pcs[it, 0]),
                 str(pcs[it, 1]),
                 str(pha[it]),
                 str(amp[it]) + '\n']
            )
            )

    # Plot RMMs
    mjp.rmm_plot(rmmfile, figname)

    # Plot composites
    filt_cube_list = [filtered_olr_cube,
                      filtered_x_wind_850_cube,
                      filtered_x_wind_200_cube,
                      filtered_precip_cube]
    composite_plots(filt_cube_list, runid, label, out_dir, rmmfile)

    # Calculate metrics
    metrics = simplified_mjo_metrics(pcs[:, 0], pcs[:, 1])
    return metrics


def composite_plots(filt_cube_list, runid, label, out_dir, rmmfile):
    ### Phase frequency plots
    out_name = runid + "_phase_freq_plot"
    figname = os.path.join(out_dir, "%s.png" % out_name)
    mjp.phasefreq_polar_plot(rmmfile, figname)

    # MJO Phase Composite Plot
    if os.path.isfile(rmmfile):
        df = pd.read_csv(rmmfile, sep=' ', names=['year', 'month', 'day', 'rmm1', 'rmm2', 'phase', 'amp'])
        # df['date'] = [datetime.datetime(y, m, d) for y, m, d in zip(df.year, df.month, df.day)]
        # Feb 29, Feb 30 is expected to throw errors as they are not real dates.
        # This is a workaround
        dates_strings = ['%s%s%s' % (str(y), str('%02d' % m), str('%02d' % d)) for y, m, d in
                         zip(df.year, df.month, df.day)]
        df['date'] = pd.DataFrame(pd.to_datetime(dates_strings, format='%Y%m%d', errors='coerce'))

    for cube in filt_cube_list:
        varname = cube.long_name

        if varname == 'x_wind':
            assert len(cube.coord('pressure').points) == 1
            pressure_level = cube.coord('pressure').points[0]
            if pressure_level == 850:
                varname = 'x_wind_850hPa'
            if pressure_level == 200:
                varname = 'x_wind_200hPa'

        # Define contour levels for plots
        if varname == 'toa_outgoing_longwave_flux':
            filt_clevels = list(range(-30, 35, 5))
            filt_colorReverse = True
        elif varname == 'precipitation_flux':
            filt_clevels = list(range(-5, 6, 1))
            filt_colorReverse = False
        elif varname == 'x_wind_850hPa':
            filt_clevels = list(range(-8, 9, 1))
            filt_colorReverse = False
        elif varname == 'x_wind_200hPa':
            filt_clevels = list(range(-14, 16, 2))
            filt_colorReverse = False

        time_coord = cube.coord('time')

        for season in ['summer', 'winter']:
            # Create results cube
            mjo_phase_composites = cube[:8]
            mjo_phase_composites.data[:, :, :] = np.nan

            out_name = '_'.join([runid, varname, season, 'phase_composites'])
            figname = os.path.join(out_dir, "%s.png" % out_name)

            fig = plt.figure(figsize=(7, 10), dpi=100)
            proj = ccrs.PlateCarree(central_longitude=180)
            cmap = getcolors('ncl_default')
            norm = colors.BoundaryNorm(filt_clevels, len(cmap.colors))

            for phase in range(1, 9):
                if season == 'summer':
                    season_df = df.loc[(df['month'].isin([5, 6, 7, 8, 9, 10]))]
                elif season == 'winter':
                    season_df = df.loc[(df['month'].isin([11, 12, 1, 2, 3, 4]))]

                phase_df = season_df.loc[(season_df['amp'] >= 1) & (season_df['phase'] == phase)]

                # get dates of the cube
                # manipulate the time string so that it allows
                # 360 day calendar which usually is not handled
                # by datetime.
                times = cube.coord('time')
                dtimes = times.units.num2date(times.points)
                dtf = pd.DataFrame(dtimes, columns=['date'])
                dtc_str = [str(dd).split()[0] for dd in dtf['date']]
                dft_cubes = pd.DataFrame(pd.to_datetime(dtc_str,
                                                        format='%Y-%m-%d', errors='coerce'),
                                         columns=['dates'])

                # Find matching indices
                phase_indices = return_indices_of_a(dft_cubes['dates'], phase_df['date'])

                ax = plt.subplot(8, 1, phase, projection=proj, facecolor='lightgrey')
                if phase_indices:
                    comp = SeasonMean_Iris(cube, phase_indices)
                    mjo_phase_composites.data[phase - 1, :, :] = comp.data

                    cf = iplt.contourf(comp, filt_clevels, cmap=cmap, norm=norm, extend='both')
                    plt.gca().coastlines()

                if phase == 1:
                    plt.title(' '.join([label, varname, season]))
                plt.text(162, 22.5, 'Phase ' + str(phase), ha='right', va='center')

            plt_ax = plt.gca()
            left, bottom, width, height = plt_ax.get_position().bounds
            first_plot_left = plt_ax.get_position().bounds[0]

            # the width of the colorbar should now be simple
            width = left - first_plot_left + width * 0.9

            # Add axes to the figure, to place the colour bar
            colorbar_axes = fig.add_axes([first_plot_left + 0.0375, bottom - 0.035, width, 0.015])

            # Add the colour bar
            cbar = plt.colorbar(cf, colorbar_axes, orientation='horizontal')
            plt.savefig(figname)
            plt.close()

            ncname = os.path.join(out_dir, "%s.nc" % out_name)
            iris.save(mjo_phase_composites, ncname)
            print('%s written.' % ncname)


def SeasonMean_Iris(var, inds):
    inds = np.array(inds)
    mean = var[inds].collapsed('time', iris.analysis.MEAN)
    return mean


def simplified_mjo_metrics(pc1, pc2):
    """Calculate Lead-lag correlation between principle components 1, and 2.
    Args:
        pc1: 1D array
        pc2: 1D array

    Returns:
        metrics dict: metric1 - 'Max corr RMM1 vs RMM2'
                      metric2 - 'Max corr RMM1 vs RMM2 lag (days)'

    Reference:
        Sperber, K. R. and Kim, D. (2012),
        Simplified metrics for the identification of the Madden - Julian oscillation
        in models.
        DOI: 10.1002/asl.378
    """
    UNDEF = -999.
    metrics = {}

    inds = np.where(pc1 != UNDEF)
    pc1 = pc1[inds]
    pc2 = pc2[inds]
    max_lag = 30
    lags = np.arange(-max_lag, max_lag + 1)
    lag_corr = mu.lagcorr(pc1, pc2, lag=lags)[:, 0]
    max_lag_corr = max(lag_corr)
    lag_of_max_corr = abs(lags[np.where(lag_corr == max_lag_corr)])[0]

    metrics['Max corr RMM1 vs RMM2'] = round(max_lag_corr, 2)
    metrics['Max corr RMM1 vs RMM2 lag (days)'] = lag_of_max_corr

    return metrics
