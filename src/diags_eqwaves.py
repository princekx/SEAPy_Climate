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

from src import data_paths
from . import mjo_plots
from src import kf_filter

def wave_amp_season(cube, out_plot_dir, runid):
    print('Dummy')
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

    print(ndjf_var.data.max())
    mjo_plots.MapPlot(ndjf_var, clevs=np.arange(0,10,1), pltname='test.png')
    plt.show()
def eqwaves_compute(cube, out_plot_dir, runid):
    print(out_plot_dir)

    cube = cube.intersection(latitude=(-30, 30), longitude=(0, 360))

    print(cube.data.max())
    # Filter transform
    ft = kf_filter.KFfilter(cube.data, spd=1)

    # MJO filter
    print('Filtering MJO...')
    cube_mjo = cube.copy()
    cube_mjo.data = ft.mjofilter()

    wave_amp_season(cube_mjo, out_plot_dir, runid)

    '''
    # Kelvin filter
    print('Filtering Kelvin...')
    cube_kelvin = cube.copy()
    cube_kelvin.data = ft.kelvinfilter()

    # Rossby wave
    print('Filtering Rossby...')
    cube_rossby = cube.copy()
    cube_rossby.data = ft.erfilter()
    '''

if __name__ == '__main__':
    eqwaves_compute(runid='obs')