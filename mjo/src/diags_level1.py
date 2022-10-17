import iris
import os
from . import mjo_utils as mu
from . import mjo_plots as mjp

def diagnos_level1(data, outplotdir, runid, label):

    varname = data.long_name
    print(varname)
    if varname == 'x_wind':
        assert len(data.coord('pressure').points) == 1
        pressure_level = data.coord('pressure').points[0]
        if pressure_level == 850:
            varname = 'x_wind_850hPa'
        if pressure_level == 200:
            varname = 'x_wind_200hPa'

    '''
    assert varname in [
            'toa_outgoing_longwave_flux',
            'precipitation_flux',
            'x_wind_850hPa',
            'x_wind_200hPa'
            ]
    '''

    print('Starting Level 1 diagnostics for %s...' % varname)

    # Define contour levels for plots
    if varname == 'toa_outgoing_longwave_flux':
        datamean_clevels = list(range(180, 280, 10))
        datamean_colorReverse = True
        datavar_clevels = list(range(200, 2800, 200))
        datavar_filt_clevels = list(range(200, 650, 50))
    elif varname == 'precipitation_flux':
        datamean_clevels = list(range(0, 16, 1))
        datamean_colorReverse = False
        datavar_clevels = list(range(0, 500, 50))
        datavar_filt_clevels = list(range(4, 44, 4))
    elif varname == 'x_wind_850hPa':
        datamean_clevels = list(range(-16, 18, 2))
        datamean_colorReverse = False
        datavar_clevels = list(range(0, 110, 10))
        datavar_filt_clevels = list(range(3, 13, 1))
    elif varname == 'x_wind_200hPa':
        datamean_clevels = list(range(-20, 44, 4))
        datamean_colorReverse = False
        datavar_clevels = list(range(0, 220, 20))
        datavar_filt_clevels = list(range(10, 70, 10))
    elif varname == 'surface_temperature':
        datamean_clevels = list(range(250, 310, 5))
        datamean_colorReverse = False
        datavar_clevels = list(range(0, 7, 1))
        datavar_filt_clevels = list(range(0, 7, 1))

    # mean plots
    summerMean = mu.SummerExtract(data)
    summerMean = summerMean.collapsed('time',iris.analysis.MEAN)

    winterMean = mu.WinterExtract(data)
    winterMean = winterMean.collapsed('time', iris.analysis.MEAN)

    ####################################################
    print('1. Plot Summer/winter mean')
    title = label + " " + varname + " mean summer"
    forename = runid + "_" + varname + "_Mean_summer"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(summerMean, title, datamean_clevels, figname, colorReverse=datamean_colorReverse)
    # save the plotted field to netcdf file
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(summerMean, ncname)


    title = label + " " + varname + " mean winter"
    forename = runid + "_" + varname + "_Mean_winter"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(winterMean, title, datamean_clevels, figname, colorReverse=datamean_colorReverse)
    # save the plotted field to netcdf file
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(winterMean, ncname)
    ####################################################


    # variance plots
    # Compute anomalies by removing a climatological annual cycle as in NCL
    harmonics = 8
    # Climatology
    clim = mu.clmDayTLL(data)
    # Smooth climatology
    clim_sm = mu.smthClmDayTLL(clim, harmonics)
    # Compute anomalies
    data = mu.calcDayAnomTLL(data, clim_sm)

    # Compute variances for both seasons
    summerVar = mu.SummerExtract(data)
    summerVar = summerVar.collapsed('time',iris.analysis.VARIANCE)

    winterVar = mu.WinterExtract(data)
    winterVar = winterVar.collapsed('time',iris.analysis.VARIANCE)

    print('2. Plot Summer/winter variance')
    title = label + " " + varname + " Var summer"
    forename = runid + "_" + varname + "_Var_summer"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(summerVar, title, datavar_clevels , figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(summerVar, ncname)

    title = label + " " + varname + " Var winter"
    forename = runid + "_" + varname + "_Var_winter"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(winterVar, title, datavar_clevels , figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(winterVar, ncname)
    ####################################################

    # filtered variance plots
    print('Filtering the series...')
    # Code below is a exact replica of the Fortran version
    # using a set of supplied filter coefficients
    datafilt = mu.Filter(data)


    print('3. Plot Summer/winter Filtered variance')
    summerFiltVar = mu.SummerExtract(datafilt)
    summerFiltVar = summerFiltVar.collapsed('time',iris.analysis.VARIANCE)

    winterFiltVar = mu.WinterExtract(datafilt)
    winterFiltVar = winterFiltVar.collapsed('time',iris.analysis.VARIANCE)

    original_varname = varname
    varname = varname + '_Filt'
    title = label + " " + varname + " Var summer"
    forename = runid + "_" + varname + "_Var_summer"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(summerFiltVar, title, datavar_filt_clevels , figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(summerFiltVar, ncname)


    title = label + " " + varname + " Var winter"
    forename = runid + "_" + varname + "_Var_winter"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(winterFiltVar, title, datavar_filt_clevels , figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(winterFiltVar, ncname)

    ####################################################
    print('4. Plot ratio of filtered variance to unfiltered variance')
    summerRatio = summerFiltVar / summerVar * 100.
    winterRatio = winterFiltVar / winterVar * 100.

    title = label + "_" + varname + " VarRatio summer"
    forename = runid + "_" + varname + "_VarRatio_summer"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(summerRatio, title, list(range(10, 60, 5)), figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(summerRatio, ncname)

    title = label + "_" + varname + " VarRatio winter"
    forename = runid + "_" + varname + "_VarRatio_winter"
    figname = os.path.join(outplotdir, "%s.png" % forename)
    mjp.MapPlot(winterRatio, title, list(range(10, 60, 5)), figname)
    ncname = os.path.join(outplotdir, "%s.nc" % forename)
    iris.save(winterRatio, ncname)
    ####################################################
    print('*' * 50)
    print(original_varname+' Level 1 diagnostics completed.')
    print('*' * 50)

    return
