import os
import sys
import numpy as np
from iris.coord_categorisation import add_year
from iris.coord_categorisation import add_month_number
from iris.coord_categorisation import add_day_of_month
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from .mycolormaps import getcolors
from src import data_paths
from src import kf_filter

def plot_wave_amps(wavename=None, varname=None, title=None, clevs=None, season=None,
                   pltname='iris_test_plot.ps', colorReverse=False, runid=None):

    # Output data/plots
    eqw_index_out_dir = os.path.join(data_paths.dirs('data_out_dir'), runid)
    sdfilename = os.path.join(eqw_index_out_dir, "%s_%s_%s_%s_sd.nc" % (runid, varname, wavename, season))
    cube = iris.load_cube(sdfilename)

    loni, lonf, lati, latf = [80, 180, -30, 30]

    cube = cube.intersection(latitude=(lati, latf), longitude=(loni, lonf))

    # Plot Lat-Lon maps
    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    #figname = os.path.join(cs_index_out_dir, "%s_%s_%s_composite.pdf" % (runid, cst, 'precip_winds850'))
    cf = iplt.contourf(cube, cmap='YlGnBu', levels=clevs, extend='both')

    plt.colorbar(cf, orientation='horizontal')
    #plt.ylim([lati, latf])
    #plt.xlim([loni, lonf])
    gl = ax.gridlines(draw_labels=True, color='white', linewidth=0.5)
    gl.xlabels_top = False
    gl.ylabels_right = False
    plt.title('%s %s amplitude %s' % (runid, wavename, season))
    plt.gca().coastlines()
    plt.savefig(pltname)
    # plt.show()


def compute_wave_amp_season(cube, runid=None, wavename=None):
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

    ndjf_var = cube[ndjf_inds].collapsed('time', iris.analysis.STD_DEV)
    jjas_var = cube[jjas_inds].collapsed('time', iris.analysis.STD_DEV)

    # Writing the standard deviation to files
    outfilename = os.path.join(eqw_index_out_dir, "%s_%s_%s_ndjf_sd.nc" % (runid, varname, wavename))
    iris.save(ndjf_var, outfilename, netcdf_format="NETCDF3_CLASSIC")

    outfilename = os.path.join(eqw_index_out_dir, "%s_%s_%s_jjas_sd.nc" % (runid, varname, wavename))
    iris.save(jjas_var, outfilename, netcdf_format="NETCDF3_CLASSIC")
    print('File written %s' %outfilename)

def eqwaves_compute(cube, out_plot_dir, runid):
    print(out_plot_dir)
    varname = cube.long_name

    cube = cube.intersection(latitude=(-30, 30), longitude=(0, 360))

    print(cube.data.max())
    # Filter transform
    ft = kf_filter.KFfilter(cube.data, spd=1)

    # MJO filter
    print('Filtering MJO...')
    wavename = 'MJO'
    cube_mjo = cube.copy()
    cube_mjo.data = ft.mjofilter()
    # Compute wave standard deviation
    compute_wave_amp_season(cube_mjo, runid=runid, wavename=wavename)
    # Plot data
    for season in ['ndjf', 'jjas']:
        pltname = os.path.join(out_plot_dir, "%s_%s_%s_%s_sd.pdf" % (runid, varname, wavename, season))
        plot_wave_amps(wavename=wavename, varname=varname, runid=runid, season=season,
                       pltname=pltname, clevs=np.arange(0, 0.5, 0.1))


    # Kelvin filter
    print('Filtering Kelvin...')
    wavename = 'Kelvin'
    cube_kelvin = cube.copy()
    cube_kelvin.data = ft.kelvinfilter()

    # Compute wave standard deviation
    compute_wave_amp_season(cube_kelvin, runid=runid, wavename=wavename)
    # Plot data
    for season in ['ndjf', 'jjas']:
        pltname = os.path.join(out_plot_dir, "%s_%s_%s_%s_sd.pdf" % (runid, varname, wavename, season))
        plot_wave_amps(wavename=wavename, varname=varname, runid=runid, season=season,
                       pltname=pltname, clevs=np.arange(0, 0.5, 0.1))

    # Rossby wave
    print('Filtering Rossby...')
    wavename = 'Rossby'
    cube_rossby = cube.copy()
    cube_rossby.data = ft.erfilter()

    # Compute wave standard deviation
    compute_wave_amp_season(cube_rossby, runid=runid, wavename=wavename)

    # Plot data
    for season in ['ndjf', 'jjas']:
        pltname = os.path.join(out_plot_dir, "%s_%s_%s_%s_sd.pdf" % (runid, varname, wavename, season))
        plot_wave_amps(wavename=wavename, varname=varname, runid=runid, season=season,
                       pltname=pltname, clevs=np.arange(0, 0.5, 0.1))

if __name__ == '__main__':
    eqwaves_compute(runid='obs')