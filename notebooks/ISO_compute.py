import numpy as np
import iris
import iris.quickplot as qplt
import iris.plot as iplt
from iris.coord_categorisation import add_year
from iris.coord_categorisation import add_month_number
from iris.coord_categorisation import add_day_of_month
import pandas as pd
import datetime

def create_dates_df(cube):
    cube = prepare_calendar(cube)
    cube_dates_df = pd.concat([pd.DataFrame(cube.coord('year').points),
                          pd.DataFrame(cube.coord('month_number').points),
                          pd.DataFrame(cube.coord('day_of_month').points)], axis=1)
    cube_dates_df.columns = ['year', 'month', 'day']
    return cube_dates_df

def precip_var_ratio_season(precip_cube, season_months = [5,6,7,8, 9, 10]):
    precip_cube_dates_df = create_dates_df(precip_cube)
    season_df = precip_cube_dates_df.loc[(precip_cube_dates_df['month'].isin(season_months))]
    seasonal_variance = precip_cube[season_df.index].collapsed('time', iris.analysis.VARIANCE)
    seasonal_filt_variance = precip_cube_filt[season_df.index].collapsed('time', iris.analysis.VARIANCE)
    return seasonal_variance, seasonal_filt_variance

def prepare_calendar(cube):
    # Setting up the dates on data
    if not cube.coords('year'):
        iris.coord_categorisation.add_year(cube, 'time', name='year')
    if not cube.coords('month_number'):
        iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')
    if not cube.coords('day_of_month'):
        iris.coord_categorisation.add_day_of_month(cube, 'time', name='day_of_month')
    return cube

data_file = '/project/MJO_GCSS/hadgem3/data/obs/GPM_PRECIP_SEA_24h_mean_2000_2020.nc'
data_file_filt = '/project/MJO_GCSS/hadgem3/data/obs/GPM_PRECIP_SEA_24h_mean_filt_2000_2020.nc'
cube = iris.load_cube(data_file_filt)
cube = cube.intersection(longitude=(90, 130), latitude=(-10,25))


'''
cube = prepare_calendar(cube)

cube_filt = lanczos_filter(cube, window=101, low=1./90., high=1./20.)

iris.save(cube_filt, data_file_filt)
print('Written %s' %data_file_filt)

'''