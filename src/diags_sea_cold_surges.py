import os
import numpy as np
from itertools import groupby
from iris.coord_categorisation import add_year
from iris.coord_categorisation import add_month_number
from iris.coord_categorisation import add_day_of_month
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import data_paths

def prepare_data(u850, v850, box_type):
    # Setting up the dates on data
    if not u850.coords('year'):
        iris.coord_categorisation.add_year(u850, 'time', name='year')
    if not u850.coords('month_number'):
        iris.coord_categorisation.add_month_number(u850, 'time', name='month_number')
    if not u850.coords('day_of_month'):
        iris.coord_categorisation.add_day_of_month(u850, 'time', name='day_of_month')

    if not v850.coords('year'):
        iris.coord_categorisation.add_year(v850, 'time', name='year')
    if not v850.coords('month_number'):
        iris.coord_categorisation.add_month_number(v850, 'time', name='month_number')
    if not v850.coords('day_of_month'):
        iris.coord_categorisation.add_day_of_month(v850, 'time', name='day_of_month')

    speed = u850.copy()
    speed.data = np.sqrt(u850.data ** 2 + v850.data ** 2)

    # box averages
    u850_ba = u850.intersection(latitude=(box_type[2], box_type[3]), longitude=(box_type[0], box_type[1]))
    u850_ba = u850_ba.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)

    v850_ba = v850.intersection(latitude=(box_type[2], box_type[3]), longitude=(box_type[0], box_type[1]))
    v850_ba = v850_ba.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)

    speed_ba = speed.intersection(latitude=(box_type[2], box_type[3]), longitude=(box_type[0], box_type[1]))
    speed_ba = speed_ba.collapsed(('latitude', 'longitude'), iris.analysis.MEAN)

    # Normalise
    #speed_ba_ts = (speed_ba.data - np.mean(speed_ba.data)) / np.std(speed_ba.data)
    speed_ba_ts = speed_ba.data
    return u850_ba, v850_ba, speed_ba_ts

def consecutive_cs_days(cs_points):
    # Counting number of consecutive cs days
    all_cs_count = np.array([len(list(g[1])) for g in groupby(cs_points) if g[0] == 1])
    cs_stats = all_cs_count[np.where(all_cs_count > 2)]
    return cs_stats

def print_dict(dictionary):
    keys = list(dictionary.keys())
    #keys.sort()
    for i in keys:
        print(i, dictionary[i])


def get_yyyymmdd(cube):
    if not cube.coords('year'):
        iris.coord_categorisation.add_year(cube, 'time', name='year')
    if not cube.coords('month_number'):
        iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')
    if not cube.coords('day_of_month'):
        iris.coord_categorisation.add_day_of_month(cube, 'time', name='day_of_month')

    years = cube.coord('year').points
    months = cube.coord('month_number').points
    days = cube.coord('day_of_month').points

    yyyymmdd = [years[i] * 10000 + months[i] * 100 + days[i] for i in range(len(years))]
    return yyyymmdd

def find_array_match_inds(long_array, short_array):
    return long_array.searchsorted(short_array)

def cold_surge_stats(u850, v850, runid='obs'):
    metrics = {}
    cs_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)

    allfile = open(os.path.join(cs_index_out_dir, "%s_NDJF_dates.txt" % runid), "w")
    csfile = open(os.path.join(cs_index_out_dir, "%s_CS_dates.txt" % runid), "w")
    msfile = open(os.path.join(cs_index_out_dir, "%s_MS_dates.txt" % runid), "w")
    esfile = open(os.path.join(cs_index_out_dir, "%s_ES_dates.txt" % runid), "w")
    cesfile = open(os.path.join(cs_index_out_dir, "%s_CES_dates.txt" % runid), "w")



    # 1. CS :Cold surge
    # Cold surge identification
    ######################################
    chang_box = [107, 115, 5, 10]  # CP Chang's 2nd domain
    chang_threshold = 9.0  # 10 # wind speed m/s
    cs_threshold = 0.75

    u850_ba, v850_ba, speed_ba_ts = prepare_data(u850, v850, chang_box)

    years = u850_ba.coord('year').points
    months = u850_ba.coord('month_number').points
    days = u850_ba.coord('day_of_month').points

    # 0. Write NDJF dates for seasonal means
    ###########################################
    allinds = [i for i in range(len(months)) if months[i] >= 11 or months[i] <= 2]

    for i in allinds:
        allfile.write(str(years[i] * 10000 + months[i] * 100 + days[i]) + '\n')
    allfile.close()
    ###########################################

    cs_points = np.zeros(len(speed_ba_ts.data))
    csinds = np.where((speed_ba_ts >= chang_threshold) & (u850_ba.data < 0.) & (v850_ba.data < 0.))[0]
    cs_points[csinds] = 1

    for i in csinds:
        csfile.write(str(years[i] * 10000 + months[i] * 100 + days[i]) + '\n')
        # print str(years[i]*10000+months[i]*100+days[i])+'\n'
    csfile.close()

    # 2. CES :Hattori box for cross equatorial surges
    ######################################
    hattori_box = [105, 115, -5, 5]
    hattori_threshold = -2.0  # m/s meridional wind

    u850_ha, v850_ha, speed_ha_ts = prepare_data(u850, v850, hattori_box)
    ces_points = np.zeros(len(speed_ba_ts.data))
    cesinds = np.where((speed_ba_ts >= cs_threshold) &
                       (u850_ba.data < 0.) &
                       (v850_ba.data < 0.) &
                       (v850_ha.data <= hattori_threshold))[0]

    ces_points[cesinds] = 1
    for i in cesinds:
        cesfile.write(str(years[i] * 10000 + months[i] * 100 + days[i]) + '\n')
        # print str(years[i]*10000+months[i]*100+days[i])+'\n'
    cesfile.close()

    # 3. MS: MERIDIONAL SURGE (COLD SURGE)
    ######################################
    # Chang et. al. (2005)
    # Averaged meridional wind between110° and 117.5 °E, along 15 °N.
    ms_box = [110, 117.5, 14, 16]
    ms_es_threshold = -8.0
    u850_ba, v850_ba, speed_ba_ts = prepare_data(u850, v850, ms_box)
    ms_points = np.zeros(len(speed_ba_ts.data))
    msinds = np.where(v850_ba.data < ms_es_threshold)[0]
    ms_points[msinds] = 1

    plt.plot(v850_ba.data)
    plt.plot(ms_points*-8)
    plt.show()

    for i in msinds:
        msfile.write(str(years[i] * 10000 + months[i] * 100 + days[i]) + '\n')
        # print str(years[i]*10000+months[i]*100+days[i])+'\n'
    msfile.close()

    #
    # 4. ES: EASTERLY SURGE
    ######################################
    # Hai et al. (2017)
    # Averaged zonal wind between 7.5° and 15 °N, along 120 °E.
    es_box = [119, 121, 7.5, 15]
    ms_es_threshold = -8.0
    u850_ba, v850_ba, speed_ba_ts = prepare_data(u850, v850, es_box)
    es_points = np.zeros(len(speed_ba_ts.data))
    esinds = np.where(u850_ba.data < ms_es_threshold)[0]
    es_points[esinds] = 1

    for i in esinds:
        esfile.write(str(years[i] * 10000 + months[i] * 100 + days[i]) + '\n')
        # print str(years[i]*10000+months[i]*100+days[i])+'\n'
    esfile.close()


    # Counting number of consecutive ces days
    cs_stats = consecutive_cs_days(cs_points)

    # Counting number of consecutive ces days
    ces_stats = consecutive_cs_days(ces_points)

    # Counting number of consecutive ms days
    ms_stats = consecutive_cs_days(ms_points)

    # Counting number of consecutive es days
    es_stats = consecutive_cs_days(es_points)

    print(len(cs_stats))
    print(len(es_stats))
    print(len(ms_stats))
    print(len(ces_stats))

    # Metric 1: average length of cold surges
    cs_mean_length = np.mean(cs_stats)
    # Metric 1: average length of CES
    ces_mean_length = np.mean(ces_stats)
    # Metric 1: average length of MS
    ms_mean_length = np.mean(ms_stats)
    # Metric 1: average length of ES
    es_mean_length = np.mean(es_stats)

    # Metric 2: average number of cold surge events in a season
    cs_mean_number = len(cs_stats) / float((years[-1] - years[0] + 1))
    ces_mean_number = len(ces_stats) / float((years[-1] - years[0] + 1))
    ms_mean_number = len(ms_stats) / float((years[-1] - years[0] + 1))
    es_mean_number = len(es_stats) / float((years[-1] - years[0] + 1))
    # print cs_mean_number #, len(cs_stats), (years[-1] - years[0]+1)

    metrics['Average length of Cold Surges (CS)'] = round(cs_mean_length, 2)
    metrics['Average length of Cross-Equatorial Surges (CES)'] = round(ces_mean_length, 2)
    metrics['Average length of Meridional Surges (MS)'] = round(ms_mean_length, 2)
    metrics['Average length of Easterly Surges (ES)'] = round(es_mean_length, 2)

    metrics['Average number of CS in a season'] = round(cs_mean_number, 2)
    metrics['Average number of CES in a season'] = round(ces_mean_number, 2)
    metrics['Average number of MS in a season'] = round(ms_mean_number, 2)
    metrics['Average number of ES in a season'] = round(es_mean_number, 2)

    print('Output written to %s' %cs_index_out_dir)

    print_dict(metrics)

    return metrics

def cold_surge_composites(vars, cstype=['NDJF', 'CS', 'CES', 'MS', 'ES'], runid='obs'):

    # Output data/plots
    cs_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)

    for cube in vars:
        for cst in cstype:
            print('Computing composite of %s' %cst)
            cs_dates_file = os.path.join(cs_index_out_dir, "%s_%s_dates.txt" % (runid, cst))
            if os.path.exists(cs_dates_file):
                allfile = open(cs_dates_file, "r")
                lines = allfile.readlines()
                cs_dates = [int(line.split()[0]) for line in lines]
                # if runid == 'obs' or (runid == 'u-ab674' and (varname == 'PRECIP' or varname == 'SST')):
                #     # ERA reanlysis has a p coordinate
                #     cube = cube[:, 0]
                # if cube.coords('t'):
                #     cube.coord('t').rename('time')
                # print
                # cube
                #
                all_dates = get_yyyymmdd(cube)
                inds = find_array_match_inds(np.array(all_dates), np.array(cs_dates))
                comp = cube[inds].collapsed('time', iris.analysis.MEAN)
                varname = cube.long_name
                outfilename = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.nc" % (runid, cst, varname))
                iris.save(comp, outfilename, netcdf_format="NETCDF3_CLASSIC")

def plot_cold_surge_composites(cstype=['CS', 'CES', 'MS', 'ES'], runid='obs'):

    # Output data/plots
    cs_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
    print(cs_index_out_dir)
    ndjf_precip_filename = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.nc" % (runid, 'NDJF', 'precipitation_flux'))
    ndjf_sst_filename = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.nc" % (runid, 'NDJF', 'surface_temperature'))
    ndjf_u850_filename = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.nc" % (runid, 'NDJF', 'x_wind_850hPa'))
    ndjf_v850_filename = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.nc" % (runid, 'NDJF', 'y_wind_850hPa'))

    if os.path.exists(ndjf_precip_filename):
        ndjf_precip_mean = iris.load_cube(ndjf_precip_filename)
    if os.path.exists(ndjf_sst_filename):
        ndjf_sst_mean = iris.load_cube(ndjf_sst_filename)
    if os.path.exists(ndjf_u850_filename):
        ndjf_u850_mean = iris.load_cube(ndjf_u850_filename)
    if os.path.exists(ndjf_v850_filename):
        ndjf_v850_mean = iris.load_cube(ndjf_v850_filename)


    for cst in cstype:
        cst_precip_filename = os.path.join(cs_index_out_dir,
                                            "%s_%s_%s_composite.nc" % (runid, cst, 'precipitation_flux'))
        cst_sst_filename = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.nc" % (runid, cst, 'surface_temperature'))
        cst_u850_filename = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.nc" % (runid, cst, 'x_wind_850hPa'))
        cst_v850_filename = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.nc" % (runid, cst, 'y_wind_850hPa'))

        if os.path.exists(cst_precip_filename):
            cst_precip_mean = iris.load_cube(cst_precip_filename)
            cst_precip_mean -= ndjf_precip_mean

        if os.path.exists(cst_sst_filename):
            cst_sst_mean = iris.load_cube(cst_sst_filename)
            cst_sst_mean -= ndjf_sst_mean

        if os.path.exists(cst_u850_filename):
            cst_u850_mean = iris.load_cube(cst_u850_filename)
            cst_u850_mean -= ndjf_u850_mean

        if os.path.exists(cst_v850_filename):
            cst_v850_mean = iris.load_cube(cst_v850_filename)
            cst_v850_mean -= ndjf_v850_mean

        x = cst_u850_mean.coord('longitude').points
        y = cst_u850_mean.coord('latitude').points
        u, v = cst_u850_mean.data, cst_v850_mean.data

        # PRECIP and winds plot
        plt.figure()
        ax = plt.axes(projection=ccrs.PlateCarree())
        figname = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.pdf" % (runid, cst, 'precip_winds850'))
        cf = iplt.contourf(cst_precip_mean, cmap='RdBu', levels=np.arange(-5, 6, 1), extend='both')

        plt.colorbar(cf, orientation='vertical')
        q = plt.quiver(x, y, u, v, units='width')
        ax.quiverkey(q, X=132, Y=-12, U=5, coordinates='data',
                     label='5 m/s', labelpos='E')
        plt.ylim([-10, 25])
        plt.xlim([95, 130])
        gl = ax.gridlines(draw_labels=True, color='white')
        gl.xlabels_top = False
        gl.ylabels_right = False
        plt.title('%s %s UV850, PRECIP' %(runid, cst))
        plt.gca().coastlines()
        plt.savefig(figname)
        #plt.show()

        # SST and winds plot
        figname = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.pdf" % (runid, cst, 'sst_winds850'))
        plt.figure()
        ax = plt.axes(projection=ccrs.PlateCarree())
        cf = iplt.contourf(cst_sst_mean, cmap='RdBu_r', levels=np.arange(-1, 1, 0.1), extend='both')

        plt.colorbar(cf, orientation='vertical')
        q = plt.quiver(x, y, u, v, units='width')
        ax.quiverkey(q, X=132, Y=-12, U=5, coordinates='data',
                     label='5 m/s', labelpos='E')
        plt.ylim([-10, 25])
        plt.xlim([95, 130])
        gl = ax.gridlines(draw_labels=True, color='white')
        gl.xlabels_top = False
        gl.ylabels_right = False
        plt.title('%s %s UV850, SST' % (runid, cst))
        plt.gca().coastlines()
        plt.savefig(figname)
        # plt.show()

if __name__ == '__main__':
    plot_cold_surge_composites(cstype=['CS', 'CES', 'MS', 'ES'], runid='obs')