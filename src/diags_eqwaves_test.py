import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib as mpl

# use Agg backend for running without X-server
mpl.use('Agg')

import iris
from iris.coord_categorisation import add_year
from iris.coord_categorisation import add_month_number
from iris.coord_categorisation import add_day_of_month
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from src import data_paths
from src import kf_filter
from scipy.signal import find_peaks

def _contour_info(varname, wavename, runid=None):
    self = dict()
    if varname == 'precipitation_flux':
        self['shade_levels'] = np.arange(0, 44, 4)
        self['contour_levels'] = np.arange(0, 0.5, 0.05)
        if wavename == 'ALL':
            self['std_levels'] = np.arange(1, 40., 5)
        elif wavename == 'MJO':
            self['std_levels'] = np.arange(1, 4., 0.25)
            self['anom_levels'] = np.linspace(-4, 4., 11)
        elif wavename == 'Kelvin':
            self['std_levels'] = np.arange(1, 4., 0.25)
            self['anom_levels'] = np.linspace(-4, 4., 11)
        else:
            self['std_levels'] = np.arange(1, 3., 0.25)

    if varname == 'x_wind_850hPa':
        self['shade_levels'] = np.arange(0, 2.2, 0.2)
        self['contour_levels'] = np.arange(0, 2.2, 0.2)
        if wavename == 'ALL':
            self['std_levels'] = np.arange(1, 5, 0.25)
        elif wavename == 'MJO':
            self['std_levels'] = np.arange(1, 3., 0.25)
            self['anom_levels'] = np.linspace(-2, 2., 11)
        elif wavename == 'Kelvin':
            self['std_levels'] = np.arange(0.5, 1.5, 0.1)
            self['anom_levels'] = np.linspace(-2, 2., 11)
        else:
            self['std_levels'] = np.arange(0.5, 1.5, 0.1)

    if varname == 'y_wind_850hPa':
        self['shade_levels'] = np.arange(0, 2.2, 0.2)
        self['contour_levels'] = np.arange(0, 2.2, 0.2)

        if wavename == 'ALL':
            self['std_levels'] = np.arange(1, 5, 0.25)
        else:
            self['std_levels'] = np.arange(0.1, 1., 0.1)
            self['anom_levels'] = np.linspace(-1, 1., 11)

    if varname == 'x_wind_200hPa':
        self['shade_levels'] = np.arange(0, 2.2, 0.2)
        self['contour_levels'] = np.arange(0, 2.2, 0.2)
        if wavename == 'ALL':
            self['std_levels'] = np.arange(1, 5, 0.25)
        elif wavename == 'MJO':
            self['std_levels'] = np.arange(1, 3., 0.25)
            self['anom_levels'] = np.linspace(-2, 2., 11)
        elif wavename == 'Kelvin':
            self['std_levels'] = np.arange(0.5, 1.5, 0.1)
            self['anom_levels'] = np.linspace(-2, 2., 11)
        else:
            self['std_levels'] = np.arange(0.5, 1.5, 0.1)
    return self


def plot_wave_std_maps(wavename=None, varname=None, title=None, season=None,
                       pltname='iris_test_plot.ps', colorReverse=False, runid=None):
    # Output data/plots
    eqw_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)

    sdfilename = os.path.join(eqw_index_out_dir, "%s_%s_%s_%s_variance.nc" % (runid, varname, wavename, season))
    if os.path.exists(sdfilename):
        cube = iris.load_cube(sdfilename)

        # Compute standard deviation
        cube.data = np.sqrt(cube.data)

    loni, lonf, lati, latf = [90, 140, -10, 20]

    cube = cube.intersection(latitude=(lati, latf), longitude=(loni, lonf))
    varname = cube.long_name
    con = _contour_info(varname, wavename)
    clevels = con['std_levels']
    print(clevels)

    # Plot Lat-Lon maps
    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    # figname = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.pdf" % (runid, cst, 'precip_winds850'))
    cf = iplt.contourf(cube, cmap='YlOrRd', extend='both', levels=clevels)  # con['contour_levels']

    plt.colorbar(cf, orientation='horizontal')
    # plt.ylim([lati, latf])
    # plt.xlim([loni, lonf])
    gl = ax.gridlines(draw_labels=True, color='white', linewidth=0.25, alpha=0.5)
    gl.xlabels_top = False
    gl.ylabels_right = False
    plt.title('%s %s %s amplitude %s' % (runid, varname, wavename, season))
    plt.gca().coastlines()
    plt.savefig(pltname, bbox_inches='tight', pad_inches=0)
    plt.close()


def plot_wave_variance_percent(wavename=None, varname=None, title=None, season=None,
                               pltname='iris_test_plot.ps', colorReverse=False, runid=None):
    # Output data/plots
    eqw_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)

    sdfilename = os.path.join(eqw_index_out_dir, "%s_%s_%s_%s_variance.nc" % (runid, varname, wavename, season))
    if os.path.exists(sdfilename):
        cube = iris.load_cube(sdfilename)

    sd_all_filename = os.path.join(eqw_index_out_dir, "%s_%s_%s_%s_variance.nc" % (runid, varname, 'ALL', season))
    try:
        cube_all = iris.load_cube(sd_all_filename)
        # finding the ratio of variances
        # avoid division by very small values
        # mask all values smaller than 0.5% of the data range
        threshold = (cube_all.data.max() - cube_all.data.min()) * 0.01

        cube_all.data = ma.masked_less(cube_all.data, threshold)

        print('Computing the ratio between %s and total variance' % wavename)
        cube.data = 100. * cube.data / cube_all.data

    except FileNotFoundError:
        print("%s not found." % sd_all_filename)
    # loni, lonf, lati, latf = [60, 180, -20, 20]
    loni, lonf, lati, latf = [90, 140, -10, 20]

    cube = cube.intersection(latitude=(lati, latf), longitude=(loni, lonf))
    varname = cube.long_name
    print(varname)
    con = _contour_info(varname, wavename)
    print(con)
    # Plot Lat-Lon maps
    clevels = [2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45]
    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    # figname = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.pdf" % (runid, cst, 'precip_winds850'))
    cf = iplt.contourf(cube, cmap='ocean_r', extend='both', levels=clevels)  # con['contour_levels']

    plt.colorbar(cf, orientation='horizontal')
    # plt.ylim([lati, latf])
    # plt.xlim([loni, lonf])
    gl = ax.gridlines(draw_labels=True, color='white', linewidth=0.25, alpha=0.5)
    gl.xlabels_top = False
    gl.ylabels_right = False
    plt.title('%s %s %s amplitude %s' % (runid, varname, wavename, season))
    plt.gca().coastlines()
    plt.savefig(pltname, bbox_inches='tight', pad_inches=0)
    plt.close()
    # plt.show()


def plot_hov_composite(wavename=None, varname=None, runid=None, title=None,
                       pltname='iris_test_plot.ps', colorReverse=False):
    # Output data/plots
    eqw_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
    sdfilename = os.path.join(eqw_index_out_dir, "%s_%s_%s_refEquator_composite.nc" % (runid, varname, wavename))

    cube = iris.load_cube(sdfilename)

    if wavename == 'Kelvin':
        cube = cube.intersection(latitude=(-5, 5), longitude=(40, 180))
    else:
        cube = cube.intersection(latitude=(-10, 10), longitude=(40, 180))

    hov_cube = cube.collapsed('latitude', iris.analysis.MEAN)
    times  = hov_cube.coord('time').points
    ntime = len(times)
    times = np.linspace(-ntime/2, ntime/2, ntime)
    #print(times)
    lons = hov_cube.coord('longitude').points

    L, T = np.meshgrid(lons, times)

    con = _contour_info(varname, wavename)
    clevels = con['anom_levels']
    cmap = 'RdBu'
    #norm = colors.BoundaryNorm(clevels, len(cmap.colors))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    CS = plt.contourf(L, T, hov_cube.data, levels=clevels,
                      cmap=cmap)
    plt.colorbar(CS, orientation='horizontal')
    #CS = plt.contour(L, T, hov_cube.data, levels=clevels, colors='k', extent='both')
    #ax.set_aspect('equal')
    # set axes range
    plt.xlim(min(lons), max(lons))
    if wavename == 'Kelvin':
        plt.ylim(min(times)/2, max(times)/2)
    else:
        plt.ylim(min(times), max(times))
    plt.title('%s %s %s [10S-10N] composite' % (runid, varname, wavename))
    plt.xlabel('Longitude (degrees east)')
    plt.ylabel('Lag/lead (days)')
    plt.grid()
    plt.savefig(pltname)#, bbox_inches='tight')
    print('Plotted %s' %pltname)
    #plt.show()
    plt.close()


def compute_wave_composite(cube, lat_range=(-10, 10),lon_range=(80, 100), wavename=None, varname=None, runid=None):
    # Output data/plots
    eqw_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
    sdfilename = os.path.join(eqw_index_out_dir, "%s_%s_%s_refEquator_composite.nc" % (runid, varname, wavename))

    ntime, nlat, nlon = cube.shape

    # Reference time series
    cube_ts = cube.intersection(latitude=lat_range, longitude=lon_range).collapsed(['latitude', 'longitude'],
                                                                                   iris.analysis.MEAN)
    # Normalise
    cube_ts_normal = (cube_ts - cube_ts.collapsed('time', iris.analysis.MEAN)) \
                     / cube_ts.collapsed('time', iris.analysis.STD_DEV)

    # Find peaks above 1 sd
    peaks, _ = find_peaks(cube_ts_normal.data, height=1.)
    npeaks = len(peaks)

    lag = 50  # days to be extracted around the peak (total window is -lag to +lag)
    comp = np.ndarray((npeaks, 2 * lag + 1, nlat, nlon), dtype=float)
    comp_cube = cube[:(2 * lag + 1)]

    for npk, peak in enumerate(peaks):
        if peak >= lag and peak <= ntime - lag - 1:
            # Append the composite
            comp[npk] = cube[peak - lag:peak + lag + 1].data

    # Composite of all events
    comp_cube.data = np.mean(comp, axis=0)

    iris.save(comp_cube, sdfilename, netcdf_format="NETCDF3_CLASSIC")
    print('Written %s .' % sdfilename)

def compute_wave_variance_season(cube, runid=None, wavename=None, metrics=None):
    if metrics == None:
        metrics = {}

    # Output data/plots
    eqw_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)

    varname = cube.long_name

    if not cube.coords('year'):
        iris.coord_categorisation.add_year(cube, 'time', name='year')
    if not cube.coords('month_number'):
        iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')
    if not cube.coords('day_of_month'):
        iris.coord_categorisation.add_day_of_month(cube, 'time', name='day_of_month')

    months = cube.coord('month_number').points

    # NDJF/JJAS dates for seasonal standard deviations
    ###########################################
    ndjf_inds = [i for i in range(len(months)) if months[i] >= 11 or months[i] <= 2]
    jjas_inds = [i for i in range(len(months)) if months[i] >= 6 and months[i] <= 9]

    print(len(ndjf_inds), len(jjas_inds))

    ndjf_var = cube[ndjf_inds].collapsed('time', iris.analysis.VARIANCE)
    jjas_var = cube[jjas_inds].collapsed('time', iris.analysis.VARIANCE)

    area_avg_variance_ndjf = ndjf_var.intersection(latitude=(-10, 12), longitude=(100, 130))
    area_avg_variance_ndjf = area_avg_variance_ndjf.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)
    metrics['Average NDJF %s %s variance' % (varname, wavename)] = round(np.asscalar(area_avg_variance_ndjf.data), 2)

    area_avg_variance_jjas = jjas_var.intersection(latitude=(-10, 12), longitude=(100, 130))
    area_avg_variance_jjas = area_avg_variance_jjas.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)
    metrics['Average JJAS %s %s variance' % (varname, wavename)] = round(np.asscalar(area_avg_variance_jjas.data), 2)

    # Writing the standard deviation to files
    outfilename = os.path.join(eqw_index_out_dir, "%s_%s_%s_ndjf_variance.nc" % (runid, varname, wavename))
    iris.save(ndjf_var, outfilename, netcdf_format="NETCDF3_CLASSIC")

    outfilename = os.path.join(eqw_index_out_dir, "%s_%s_%s_jjas_variance.nc" % (runid, varname, wavename))
    iris.save(jjas_var, outfilename, netcdf_format="NETCDF3_CLASSIC")
    print('File written %s' % outfilename)

    print(metrics)
    return metrics


def eqwaves_compute(cube, out_plot_dir, runid, compute=True, write_wave_data=True,
                    write_wave_data_dir=None):
    metrics = {}
    print(out_plot_dir)
    varname = cube.long_name

    cube = cube.intersection(latitude=(-30, 30), longitude=(0, 360))

    print(cube.data.max())
    if not compute:
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('NO COMPUTATIONS WILL BE PERFORMED. ONLY PLOTTING...')
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

    # Compute seasonal variance first
    wavename = 'ALL'
    if compute:
        cube_all = cube.copy()
        # Compute wave standard deviation
        metrics.update(compute_wave_variance_season(cube_all, runid=runid, wavename=wavename, metrics=metrics))

    # Plot standard deviation maps
    for season in ['ndjf', 'jjas']:
        pltname = os.path.join(out_plot_dir, "%s_%s_%s_%s_stddev.pdf" % (runid, varname, wavename, season))
        plot_wave_std_maps(wavename=wavename, varname=varname, runid=runid, season=season,
                           pltname=pltname)

        print('Plotted %s' % pltname)

    # MJO filter
    print('Filtering MJO...')
    wavename = 'MJO'
    if compute:
        # Filter transform
        ft = kf_filter.KFfilter(cube.data, spd=1)

        cube_mjo = cube.copy()
        cube_mjo.data = ft.mjofilter()

        print('MJO Filtered')

        # Compute equatorial composite
        try:
            compute_wave_composite(cube_mjo, lat_range=(-10, 10),lon_range=(80, 100),
                                   wavename=wavename, varname=varname, runid=runid)
            # Plot hovmoller
            plt_hov_name = os.path.join(out_plot_dir, "%s_%s_%s_composite_hovmoller.pdf" %
                                   (runid, varname, wavename))

            plot_hov_composite(wavename=wavename, varname=varname, runid=runid, pltname=plt_hov_name)
            print('Plotted %s' % plt_hov_name)
        except:
            print('Failed in composite computation.')

        # Compute wave variance
        metrics.update(compute_wave_variance_season(cube_mjo, runid=runid,
                                                    wavename=wavename, metrics=metrics))

        # Write the wave filtered data
        if write_wave_data and write_wave_data_dir:
            wave_data_file_name = os.path.join(write_wave_data_dir,
                                               "%s_%s_%s_filtered.nc"
                                               % (runid, varname, wavename))
            iris.save(cube_mjo, wave_data_file_name, netcdf_format="NETCDF3_CLASSIC")
            print('Written %s .' % wave_data_file_name)

    # Plot data
    for season in ['ndjf', 'jjas']:
        # Plot standard deviation maps
        pltname = os.path.join(out_plot_dir, "%s_%s_%s_%s_stddev.pdf" % (runid, varname, wavename, season))
        plot_wave_std_maps(wavename=wavename, varname=varname, runid=runid, season=season,
                           pltname=pltname)
        print('Plotted %s' % pltname)

        # Plot ratio of variances
        pltname = os.path.join(out_plot_dir, "%s_%s_%s_%s_variance_percent.pdf" % (runid, varname, wavename, season))
        plot_wave_variance_percent(wavename=wavename, varname=varname, runid=runid, season=season,
                                   pltname=pltname)
        print('Plotted %s' % pltname)

    # Kelvin filter
    print('Filtering Kelvin...')
    wavename = 'Kelvin'
    if compute:
        cube_kelvin = cube.copy()
        cube_kelvin.data = ft.kelvinfilter()

        # Compute wave variance
        metrics.update(compute_wave_variance_season(cube_kelvin, runid=runid, wavename=wavename, metrics=metrics))

        try:
            compute_wave_composite(cube_kelvin, lat_range=(-10, 10), lon_range=(80, 100),
                                   wavename=wavename, varname=varname, runid=runid)
            # Plot hovmoller
            plt_hov_name = os.path.join(out_plot_dir, "%s_%s_%s_composite_hovmoller.pdf" %
                                   (runid, varname, wavename))
            plot_hov_composite(wavename=wavename, varname=varname, runid=runid, pltname=plt_hov_name)
            print('Plotted %s' % pltname)
        except:
            print('Failed in composite computation.')

    # Plot data
    for season in ['ndjf', 'jjas']:
        # Plot standard deviation maps
        pltname = os.path.join(out_plot_dir, "%s_%s_%s_%s_stddev.pdf" % (runid, varname, wavename, season))
        plot_wave_std_maps(wavename=wavename, varname=varname, runid=runid, season=season,
                           pltname=pltname)
        print('Plotted %s' % pltname)

        # Plot ratio of variances
        pltname = os.path.join(out_plot_dir, "%s_%s_%s_%s_variance_percent.pdf" % (runid, varname, wavename, season))
        plot_wave_variance_percent(wavename=wavename, varname=varname, runid=runid, season=season,
                                   pltname=pltname)
        print('Plotted %s' % pltname)

    # Rossby wave
    print('Filtering Rossby...')
    wavename = 'Rossby'
    if compute:
        cube_rossby = cube.copy()
        cube_rossby.data = ft.erfilter()

        # Compute wave variance
        metrics.update(compute_wave_variance_season(cube_rossby, runid=runid, wavename=wavename, metrics=metrics))

        # Write the wave filtered data
        if write_wave_data and write_wave_data_dir:
            wave_data_file_name = os.path.join(write_wave_data_dir,
                                               "%s_%s_%s_filtered.nc"
                                               % (runid, varname, wavename))
            iris.save(cube_rossby, wave_data_file_name, netcdf_format="NETCDF3_CLASSIC")
            print('Written %s .' % wave_data_file_name)

    # Plot data
    for season in ['ndjf', 'jjas']:
        # Plot standard deviation maps
        pltname = os.path.join(out_plot_dir, "%s_%s_%s_%s_stddev.pdf" % (runid, varname, wavename, season))
        plot_wave_std_maps(wavename=wavename, varname=varname, runid=runid, season=season,
                           pltname=pltname)
        print('Plotted %s' % pltname)

        # Plot ratio of variances
        pltname = os.path.join(out_plot_dir, "%s_%s_%s_%s_variance_percent.pdf" % (runid, varname, wavename, season))
        plot_wave_variance_percent(wavename=wavename, varname=varname, runid=runid, season=season,
                                   pltname=pltname)
        print('Plotted %s' % pltname)

    print(metrics)
    return metrics


if __name__ == '__main__':
    # data
    obs = {}
    obs['runid'] = 'obs'
    obs['start_date'] = '1989/01/01'
    obs['end_date'] = '1998/12/01'
    obs['data_retrieve_dir'] = '/project/MJO_GCSS/hadgem3/data/obs/SEAPy_data'

    run = obs
    data_root = os.path.join(run['data_retrieve_dir'], 'obs')
    out_plot_dir = '/project/MJO_GCSS/hadgem3/data/SEAPy_output'
    print(data_root)
    data_file = os.path.join(data_root, 'obs_PRECIP.pp.nc')
    cube = iris.load_cube(data_file)
    cube.long_name = 'precipitation_flux'
    #print(cube)

    wave_metrics = eqwaves_compute(cube, out_plot_dir, 'obs', compute=True,
                                   write_wave_data=True, write_wave_data_dir=data_root)
    # print(wave_metrics)
    #plot_hov_composite(wavename='MJO', varname='precipitation_flux', runid='obs')