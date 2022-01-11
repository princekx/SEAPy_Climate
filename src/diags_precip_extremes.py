import os
import datetime

import numpy as np
import matplotlib as mpl

# use Agg backend for running without X-server
mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import iris
import iris.plot as iplt
from iris.coord_categorisation import add_year
from iris.coord_categorisation import add_month_number
from iris.coord_categorisation import add_day_of_month
import cartopy.crs as ccrs
import pandas as pd
from src import data_paths


def prepare_calendar(cube):
    '''
    Function to add year, month and day categorisation coordinates
    to the cube.
    :param cube: iris cube with time coordinate
    :type cube: iris cube
    :return:
    :rtype:
    '''

    # Setting up the dates on data
    if not cube.coords('year'):
        iris.coord_categorisation.add_year(cube, 'time', name='year')
    if not cube.coords('month_number'):
        iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')
    if not cube.coords('day_of_month'):
        iris.coord_categorisation.add_day_of_month(cube, 'time', name='day_of_month')
    return cube


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


def compute_percentiles(cube, runid='obs', cstype='NDJF', percent=[95]):
    '''
    Function to compute percentile maps along time coordinate
    for a selection of season or dates aaasssss given by cs_dates_file.
    :param cube:
    :type cube:
    :param runid:
    :type runid:
    :param cstype:
    :type cstype:
    :param percent:
    :type percent:
    :return:
    :rtype:
    '''
    # Output  data / plots
    cs_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
    # read the CS/seasonal dates
    cs_dates_file = os.path.join(cs_index_out_dir, "%s_%s_dates.txt" % (runid, cstype))
    if os.path.exists(cs_dates_file):
        print(cs_dates_file)
        df = pd.read_csv(cs_dates_file, header=None)
        df.columns = ['dates']
        df_times_cst = [d for d in df.dates]

        # Dataframe method does not support 360 day calendars
        # So switiching it back to made up dates
        #df['datetime'] = pd.to_datetime(df['dates'], format='%Y%m%d')
        #dft_cst = pd.DataFrame(df['datetime'])

        # get dates of the cube and take care of 360 day calendars
        times = cube.coord('time')
        dtimes_cube = [int(time.strftime('%Y%m%d')) for time in times.units.num2date(times.points)]

        # Dataframe method does not support 360 day calendars
        # So switiching it back to made up dates
        #dft_cubes = pd.DataFrame(pd.to_datetime(dtimes, format='%Y%m%d'))
        #dft_cubes.columns = ['dates']

        # get indices of cube that matches the CS dates
        indices = return_indices_of_a(dtimes_cube, df_times_cst)
        print(len(indices))

        # compute percentiles
        percentiles = cube[indices].collapsed('time', iris.analysis.PERCENTILE, percent=percent)
        return percentiles
    else:
        print('%s not found. Skipping computation' % cs_dates_file)


def precip_extremes_cs(precip_cubes, runid='obs', percentile=95,
                       compute=True, plot=True, frac_plot=True):
    '''
    Computes, writes and/or plots the percentile maps
    :param precip_cubes:
    :type precip_cubes:
    :param runid:
    :type runid:
    :param percentile:
    :type percentile:
    :param compute:
    :type compute:
    :param plot:
    :type plot:
    :return:
    :rtype:
    '''

    cs_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)

    if compute:
        # subset region
        precip_cubes = precip_cubes.intersection(longitude=(95, 130), latitude=(-10, 25))
        precip_cubes = prepare_calendar(precip_cubes)

        for cst in ['NDJF', 'CS', 'CES']:
            precip_perc = compute_percentiles(precip_cubes, cstype=cst, percent=[percentile])
            print(precip_perc)
            outfilename = os.path.join(cs_index_out_dir, "%s_%s_%s_%s_percentile.nc" %
                                       (runid, cst, 'precip', percentile))
            iris.save(precip_perc, outfilename, netcdf_format="NETCDF3_CLASSIC")
            print('%s written.' % outfilename)

    if plot:
        for cst in ['NDJF', 'CS', 'CES']:
            filename = os.path.join(cs_index_out_dir, "%s_%s_%s_%s_percentile.nc" %
                                    (runid, cst, 'precip', percentile))
            precip_perc = iris.load_cube(filename)

            plt.figure()
            ax = plt.axes(projection=ccrs.PlateCarree())
            figname = os.path.join(cs_index_out_dir, '%s_%s_%s_%s_percentile.pdf' %
                                   (runid, cst, 'precip', percentile))

            levels = np.arange(10, 80, 10)
            cmap = plt.get_cmap('GnBu')
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            cf = iplt.pcolormesh(precip_perc, cmap=cmap, norm=norm)

            plt.colorbar(cf, orientation='vertical')
            plt.ylim([-10, 25])
            plt.xlim([95, 130])
            gl = ax.gridlines(draw_labels=True, alpha=0.5)
            gl.xlabels_top = False
            gl.ylabels_right = False
            plt.title('%s %s %s percentile PRECIP' % (runid, cst, percentile))
            plt.gca().coastlines()
            plt.savefig(figname, bbox_inches='tight', pad_inches=0)
            print('%s plotted.' % figname)

    if frac_plot:
        # seasonal values to compute fraction with
        filename = os.path.join(cs_index_out_dir, "%s_%s_%s_%s_percentile.nc" %
                                (runid, 'NDJF', 'precip', percentile))
        seasonal_precip_perc = iris.load_cube(filename)

        for cst in ['CS', 'CES']:
            filename = os.path.join(cs_index_out_dir, "%s_%s_%s_%s_percentile.nc" %
                                    (runid, cst, 'precip', percentile))
            precip_perc = iris.load_cube(filename)

            plt.figure()
            ax = plt.axes(projection=ccrs.PlateCarree())
            figname = os.path.join(cs_index_out_dir, '%s_%s_%s_%s_percentile_fraction.pdf' %
                                   (runid, cst, 'precip', percentile))
            levels = np.arange(-100, 120, 20)
            cmap = plt.get_cmap('RdBu')
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
            cf = iplt.pcolormesh(100. * (precip_perc - seasonal_precip_perc) / seasonal_precip_perc, cmap=cmap,
                                 norm=norm)
            plt.colorbar(cf, orientation='vertical')
            plt.ylim([-10, 25])
            plt.xlim([95, 130])
            gl = ax.gridlines(draw_labels=True, alpha=0.5)
            gl.xlabels_top = False
            gl.ylabels_right = False
            plt.title('%s perc change in %s percentile PRECIP \n due to %s' % (runid, percentile, cst))
            plt.gca().coastlines()
            plt.savefig(figname, bbox_inches='tight', pad_inches=0)
            print('%s plotted.' % figname)


def mjo_phase_indices(cube, runid=None, season='NDJF', phases=[1]):
    '''
    Finds the cube indices that match MJO phase input
    also returns the NDJF indices for comparison
    :param cube:
    :type cube:
    :param runid:
    :type runid:
    :param season:
    :type season:
    :param phase:
    :type phase:
    :return: season_indices, phase_indices
    :rtype:
    '''
    # mjo_dates_file = '/project/MJO_GCSS/hadgem3/data/MJOPy_output/obs/RMMs_obs.txt'
    # Output  data / plots
    mjo_index_out_dir = os.path.join(data_paths.dirs('mjo_data_out_dir'), runid)
    # read the CS/seasonal dates
    mjo_dates_file = os.path.join(mjo_index_out_dir, "RMMs_%s.txt" % runid)

    df = pd.read_csv(mjo_dates_file, sep=' ', names=['year', 'month', 'day', 'rmm1', 'rmm2', 'phase', 'amp'])
    df['date'] = [datetime.datetime(y, m, d) for y, m, d in zip(df.year, df.month, df.day)]
    if season == 'NDJF':
        ndjf_df = df.loc[(df['month'].isin([11, 12, 1, 2]))]
        phase_df = ndjf_df.loc[(ndjf_df['amp'] >= 1) & (ndjf_df['phase'].isin(phases))]
        # get dates of the cube
        times = cube.coord('time')
        dtimes = times.units.num2date(times.points)
        dft_cubes = pd.DataFrame(pd.to_datetime(dtimes, format='%Y%m%d'))
        dft_cubes.columns = ['dates']

        # get indices of cube that matches the CS dates
        season_indices = return_indices_of_a(dft_cubes['dates'], ndjf_df['date'])
        phase_indices = return_indices_of_a(dft_cubes['dates'], phase_df['date'])
        return season_indices, phase_indices


def precip_extremes_mjo(cube, runid='obs', percentile=95,
                        compute=True, plot=True, frac_plot=True):
    '''
    Function computes the percentile maps of seasonal and MJO phases
    writes out a cube of size = 9, the first element being the seasonal percentile
    and subsequent values are for the respective MJO phases
    :param cube:
    :type cube:
    :param runid:
    :type runid:
    :param percentile:
    :type percentile:
    :param compute:
    :type compute:
    :param plot:
    :type plot:
    :return:
    :rtype:
    '''

    #mjo_index_out_dir = os.path.join(data_paths.dirs('mjo_data_out_dir'), runid)
    mjo_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)

    cube = cube.intersection(longitude=(95, 130), latitude=(-10, 25))
    if compute:
        # Computing phases 2,3,4
        ndjf_indices, p234_indices = mjo_phase_indices(cube, phases=[2, 3, 4],
                                                       runid=runid, season='NDJF')
        # first cube is seasonal extremes
        seas_perc = cube[ndjf_indices].collapsed('time',
                                                 iris.analysis.PERCENTILE,
                                                 percent=[percentile])
        outfilename = os.path.join(mjo_index_out_dir, "%s_%s_%s_%s_percentile.nc" %
                                   (runid, 'NDJF', 'precip', percentile))
        iris.save(seas_perc, outfilename, netcdf_format="NETCDF3_CLASSIC")
        print('%s written.' % outfilename)

        # extremes for phases 2,3,4
        mjo_perc_234 = cube[p234_indices].collapsed('time',
                                                    iris.analysis.PERCENTILE,
                                                    percent=[percentile])
        outfilename = os.path.join(mjo_index_out_dir, "%s_%s_%s_%s_percentile.nc" %
                                   (runid, 'MJO_234', 'precip', percentile))
        iris.save(mjo_perc_234, outfilename, netcdf_format="NETCDF3_CLASSIC")
        print('%s written.' % outfilename)

        # Computing extremes for phases 6,7,8
        ndjf_indices, p678_indices = mjo_phase_indices(cube, phases=[6, 7, 8],
                                                       runid=runid, season='NDJF')
        mjo_perc_678 = cube[p678_indices].collapsed('time',
                                                    iris.analysis.PERCENTILE,
                                                    percent=[percentile])
        outfilename = os.path.join(mjo_index_out_dir, "%s_%s_%s_%s_percentile.nc" %
                                   (runid, 'MJO_678', 'precip', percentile))
        iris.save(mjo_perc_678, outfilename, netcdf_format="NETCDF3_CLASSIC")
        print('%s written.' % outfilename)

    if plot:
        for mjo_phase in ['NDJF', 'MJO_234', 'MJO_678']:
            outfilename = os.path.join(mjo_index_out_dir, "%s_%s_%s_%s_percentile.nc" %
                                       (runid, mjo_phase, 'precip', percentile))
            if os.path.exists(outfilename):
                cube = iris.load_cube(outfilename)
                plt.figure()
                ax = plt.axes(projection=ccrs.PlateCarree())
                figname = os.path.join(mjo_index_out_dir, '%s_%s_%s_%s_percentile.pdf' %
                                       (runid, mjo_phase, 'precip', percentile))

                levels = np.arange(10, 80, 10)
                cmap = plt.get_cmap('GnBu')
                norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                cf = iplt.pcolormesh(cube, cmap=cmap, norm=norm)

                plt.colorbar(cf, orientation='vertical')
                plt.ylim([-10, 25])
                plt.xlim([95, 130])
                gl = ax.gridlines(draw_labels=True, alpha=0.5)
                gl.xlabels_top = False
                gl.ylabels_right = False
                plt.title('%s %s %s percentile PRECIP' % (runid, mjo_phase, percentile))
                plt.gca().coastlines()
                plt.savefig(figname, bbox_inches='tight', pad_inches=0)
                print('%s plotted.' % figname)
    if frac_plot:
        # seasonal values to compute fraction with
        outfilename = os.path.join(mjo_index_out_dir, "%s_%s_%s_%s_percentile.nc" %
                                   (runid, 'NDJF', 'precip', percentile))
        seasonal_precip_perc = iris.load_cube(outfilename)
        for mjo_phase in ['MJO_234', 'MJO_678']:
            outfilename = os.path.join(mjo_index_out_dir, "%s_%s_%s_%s_percentile.nc" %
                                       (runid, mjo_phase, 'precip', percentile))
            if os.path.exists(outfilename):
                cube = iris.load_cube(outfilename)
                plt.figure()
                ax = plt.axes(projection=ccrs.PlateCarree())
                figname = os.path.join(mjo_index_out_dir, '%s_%s_%s_%s_percentile.pdf' %
                                       (runid, mjo_phase, 'precip', percentile))
                # cf = iplt.contourf(100. * (cube - seasonal_precip_perc) / seasonal_precip_perc,
                #                   levels=np.arange(-100, 120, 20), cmap='RdBu', extend='both')
                levels = np.arange(-100, 120, 20)
                cmap = plt.get_cmap('RdBu')
                norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
                cf = iplt.pcolormesh(100. * (cube - seasonal_precip_perc) / seasonal_precip_perc, cmap=cmap,
                                     norm=norm)

                plt.colorbar(cf, orientation='vertical')
                plt.ylim([-10, 25])
                plt.xlim([95, 130])
                gl = ax.gridlines(draw_labels=True, alpha=0.5)
                gl.xlabels_top = False
                gl.ylabels_right = False
                plt.title('%s perc change in %s percentile PRECIP \n due to %s' % (runid, percentile, mjo_phase))
                plt.gca().coastlines()
                plt.savefig(figname, bbox_inches='tight', pad_inches=0)
                print('%s plotted.' % figname)


if __name__ == '__main__':
    precip_cubes = iris.load_cube('/project/MJO_GCSS/hadgem3/data/obs/SEAPy_data/obs/obs_PRECIP.pp.nc')
    precip_extremes_cs(precip_cubes, runid='obs', percentile=95, compute=True, plot=True, frac_plot=True)
    precip_extremes_mjo(precip_cubes, runid='obs', percentile=95, compute=True, plot=True, frac_plot=True)
