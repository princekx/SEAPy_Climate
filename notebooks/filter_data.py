import numpy as np
import iris
import iris.quickplot as qplt
import iris.plot as iplt
from iris.coord_categorisation import add_year
from iris.coord_categorisation import add_month_number
from iris.coord_categorisation import add_day_of_month

def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    """
    order = ((window - 1) // 2 ) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    return w[1:-1]

def lanczos_filter(cube, window=101, low=1./100., high=1./20.):
    wgts_high = low_pass_weights(window, high)
    wgts_low = low_pass_weights(window, low)
    cube_filt = cube.copy()
    cube_filt.data = np.apply_along_axis(lambda m: (np.convolve(m, wgts_high, mode='same') -
                            np.convolve(m, wgts_low, mode='same')), axis=0, arr=cube.data)
    return cube_filt

def prepare_calendar(cube):
    # Setting up the dates on data
    if not cube.coords('year'):
        iris.coord_categorisation.add_year(cube, 'time', name='year')
    if not cube.coords('month_number'):
        iris.coord_categorisation.add_month_number(cube, 'time', name='month_number')
    if not cube.coords('day_of_month'):
        iris.coord_categorisation.add_day_of_month(cube, 'time', name='day_of_month')
    return cube


data_file = '/project/MJO_GCSS/hadgem3/data/obs/ERA5_tropical_U850_2000_2020_dailymean.nc'
data_file_filt = '/project/MJO_GCSS/hadgem3/data/obs/ERA5_U850_SEA_24h_mean_2000_2020.nc'
cube = iris.load_cube(data_file)
cube = cube.intersection(longitude=(90, 130), latitude=(-10,25))

cube = prepare_calendar(cube)

#cube_filt = lanczos_filter(cube, window=101, low=1./90., high=1./20.)

#iris.save(cube_filt, data_file_filt)
iris.save(cube, data_file_filt)
print('Written %s' %data_file_filt)

