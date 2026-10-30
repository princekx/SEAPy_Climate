import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from .plot_utils import save_figure
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import matplotlib.colors as colors
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Circle
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import os, sys
import pandas as pd
from .mycolormaps import getcolors
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import textwrap

try:  # seaborn only styles the plots, so it stays optional
    import seaborn as sns
except ImportError:
    sns = None


def _set_plot_theme():
    if sns is not None:
        sns.set_theme(style="ticks")


#mpl.rc('font', family='Times New Roman')
_set_plot_theme()

def MapPlot(var2plot, title=None, clevs=None, pltname='iris_test_plot.ps', colorReverse=False):
    # Plot Lat-Lon maps
    fig = plt.figure(figsize=(8, 3), dpi=100)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    if colorReverse:
        cf = iplt.contourf(var2plot, clevs, cmap=getcolors('WhBlGrYeRe', colorReverse=True), extend='both')
    else:
        cf = iplt.contourf(var2plot, clevs, cmap=getcolors('WhBlGrYeRe'), extend='both')

    #iplt.pcolormesh(var2plot, cmap=getcolors('WhBlGrYeRe') )
    ax.coastlines()
    ax.gridlines(alpha=0.0)
    #ax.set_global()
    if len(title) > 40:
        title = textwrap.fill(title, width=40, break_long_words=False,
                                break_on_hyphens=False)
    plt.title(title)

    ax.set_xticks([0, 60, 120, 180, 240, 300, 360], crs=ccrs.PlateCarree())
    ax.set_yticks([-30, 0, 30], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=0, fontsize=10)
    labels = ax.get_yticklabels()
    plt.setp(labels, rotation=0, fontsize=10)

    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    first_plot_left = plt_ax.get_position().bounds[0]

    # the width of the colorbar should now be simple
    width = left - first_plot_left + width * 0.9

    # Add axes to the figure, to place the colour bar
    #colorbar_axes = fig.add_axes([first_plot_left + 0.0375, bottom + 0.1, width, 0.05])
    colorbar_axes = fig.add_axes([first_plot_left + 0.0375, bottom - 0.2, width, 0.05])
    # Add the colour bar
    cbar = plt.colorbar(cf, colorbar_axes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=10)

    save_figure(pltname, bbox_inches='tight')
    #plt.show()
    plt.close()

def HovTimeLon(var, time, lons, levels=None, title=None, figname=None):
    L, T = np.meshgrid(lons, time)
    if len(title) > 40:
        title = textwrap.fill(title, width=40, break_long_words=False,
                                break_on_hyphens=False)
    cmap = getcolors('ncl_default')
    norm = colors.BoundaryNorm(levels, len(cmap.colors))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    CS = plt.contourf(L, T, var, levels=levels,
                      cmap=cmap, norm=norm)
    plt.colorbar(CS, orientation='horizontal')
    CS = plt.contour(L, T, var, levels=levels, colors='k')
    #ax.set_aspect('equal')
    # set axes range
    plt.xlim(min(lons), max(lons))
    plt.ylim(min(time), max(time))
    plt.title(title)#
    plt.xlabel('Longitude (degrees east)')
    plt.ylabel('Lag/lead (days)')
    plt.grid()
    save_figure(figname, bbox_inches='tight')
    #plt.show()
    plt.close()

def PlotSpecRaw(spec, freq, wave, levels=np.linspace(-2, 2, 21), \
             title='', figname='spec_test.ps'):

    if len(title) > 40:
        title = textwrap.fill(title, width=40, break_long_words=False,
                                break_on_hyphens=False)
    plt.clf()
    minfrq4plt = 0.
    maxfrq4plt = 0.8
    minwav4plt = -15
    maxwav4plt = 15

    minfrq = minfrq4plt
    maxfrq = min([maxfrq4plt, max(freq)])
    F, W = np.meshgrid(wave, freq)

    cmap = getcolors('ncl_default')
    norm = colors.BoundaryNorm(levels, len(cmap.colors))
    norm = colors.Normalize(clip=False)
    CS = plt.contourf(F, W, spec, levels=levels,
                      cmap=cmap, norm=norm, extend='both')
    bar = plt.colorbar(CS)

    tick_locs = levels
    tick_labels = levels
    bar.locator = ticker.FixedLocator(tick_locs)
    bar.formatter = ticker.FixedFormatter(tick_labels)
    bar.update_ticks()

    # set axes range
    plt.xlim(minwav4plt, maxwav4plt)
    plt.ylim(minfrq, maxfrq)

    # Lines
    plt.plot([0, 0], [0, 0.5], 'k--', lw=0.5)

    # Line markers of periods
    frqs = [80, 30, 6, 3]
    for frq in frqs:
        plt.plot([-15, 15], [1. / frq, 1. / frq], 'k--', lw=0.5)
        plt.text(-14.7, 1. / frq, str(frq) + ' days', {'color' : 'k'})

    plt.title(title)#
    plt.xlabel('Westward     Zonal Wave Number     Eastward')
    plt.ylabel('Frequency (CPD)')

    save_figure(figname)#,bbox_inches='tight')
    #plt.show()
    plt.close()

def plotAntiSymmetric(spec, freq, wave, Apzwn, \
         Afreq, levels=np.linspace(-2, 2, 21), \
         title='', figname='specAntiSym_test.ps'):

    if len(title) > 40:
        title = textwrap.fill(title, width=40, break_long_words=False,
                                break_on_hyphens=False)
    plt.clf()
    minfrq4plt = 0.
    maxfrq4plt = 0.8
    minwav4plt = -15
    maxwav4plt = 15

    minfrq = minfrq4plt
    maxfrq = min([maxfrq4plt, max(freq)])
    F, W = np.meshgrid(wave, freq)

    cmap = getcolors('ncl_default')
    norm = colors.BoundaryNorm(levels, len(cmap.colors))

    CS = plt.contourf(F, W, spec, levels=levels,
                      cmap=cmap, norm=norm, extend='both')
    bar = plt.colorbar(CS)

    tick_locs = levels
    tick_labels = levels
    bar.locator = ticker.FixedLocator(tick_locs)
    bar.formatter = ticker.FixedFormatter(tick_labels)
    bar.update_ticks()

    # set axes range
    plt.xlim(minwav4plt, maxwav4plt)
    plt.ylim(minfrq, maxfrq)

    # Lines
    plt.plot([0, 0], [0, 0.5], 'k--', lw=0.5)

    # Line markers of periods
    frqs = [80, 30, 6, 3]
    for frq in frqs:
        plt.plot([-15, 15], [1. / frq, 1. / frq], 'k--', lw=0.5)
        plt.text(-14.7, 1. / frq, str(frq) + ' days', {'color' : 'k'})

    plt.title(title)#
    plt.xlabel('Westward     Zonal Wave Number     Eastward')
    plt.ylabel('Frequency (CPD)')

    # Symmetric waves
    # Equatorial Rossby
    for i in range(3):
        for j in range(3):
            plt.plot(Apzwn[i, j, :], Afreq[i, j, :], 'k', lw=0.5)

    plt.text(-10., 0.15, "MRG", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(-3., 0.58, "n=2 IG", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(6., 0.4, "n=0 EIG", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(-3., 0.475, "h=12", {'color' : 'k', 'backgroundcolor':'w'})

    save_figure(figname)#,bbox_inches='tight')
    #plt.show()
    plt.close()

def plotSymmetric(spec, freq, wave, Apzwn, \
         Afreq, levels=np.linspace(-2, 2, 21), \
         title='', figname='specSym_test.ps'):

    if len(title) > 40:
        title = textwrap.fill(title, width=40, break_long_words=False,
                                break_on_hyphens=False)
    plt.clf()
    minfrq4plt = 0.
    maxfrq4plt = 0.8
    minwav4plt = -15
    maxwav4plt = 15

    minfrq = minfrq4plt
    maxfrq = min([maxfrq4plt, max(freq)])
    F, W = np.meshgrid(wave, freq)

    cmap = getcolors('ncl_default')
    norm = colors.BoundaryNorm(levels, len(cmap.colors))

    CS = plt.contourf(F, W, spec, levels=levels,
                      cmap=cmap, norm=norm, extend='both')
    bar = plt.colorbar(CS)

    tick_locs = levels
    tick_labels = levels
    bar.locator = ticker.FixedLocator(tick_locs)
    bar.formatter = ticker.FixedFormatter(tick_labels)
    bar.update_ticks()

    # set axes range
    plt.xlim(minwav4plt, maxwav4plt)
    plt.ylim(minfrq, maxfrq)

    # Lines
    plt.plot([0, 0], [0, 0.5], 'k--', lw=0.5)

    # Line markers of periods
    frqs = [80, 30, 6, 3]
    for frq in frqs:
        plt.plot([-15, 15], [1. / frq, 1. / frq], 'k--', lw=0.5)
        plt.text(-14.7, 1. / frq, str(frq) + ' days', {'color' : 'k'})

    plt.title(title)#
    plt.xlabel('Westward     Zonal Wave Number     Eastward')
    plt.ylabel('Frequency (CPD)')

    # Symmetric waves
    # Equatorial Rossby
    for i in range(3, 6):
        for j in range(3):
            plt.plot(Apzwn[i, j, :], Afreq[i, j, :], 'k', lw=0.5)

    plt.text(11.5, 0.4, "Kelvin", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(-10.7, 0.07, "n=1 ER", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(-3., 0.45, "n=1 IG", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(-14., 0.46, "h=12", {'color' : 'k', 'backgroundcolor':'w'})

    save_figure(figname)#,bbox_inches='tight')
    #plt.show()
    plt.close()

def mjo_wavenum_freq_season_plot(pow_cube, levels=np.arange(0, 3, 0.2), \
             title='', figname='wavenum_freq_season_plot.ps'):

    if len(title) > 40:
        title = textwrap.fill(title, width=40, break_long_words=False,
                                break_on_hyphens=False)
    NW = 6
    fMin = -0.05
    fMax = 0.05

    constraint = iris.Constraint(frequency=lambda cell: fMin <= cell <= fMax, \
                                 wavenumber=lambda cell: 0 <= cell <= NW)
    pow_cube = pow_cube.extract(constraint)

    freq = pow_cube.coord('frequency').points
    wave = pow_cube.coord('wavenumber').points

    W, F = np.meshgrid(freq, wave)

    # set zeroth frequency to minimum value
    izero = np.where(freq == 0)[0][0]
    pow_cube.data[:,izero] = np.min(pow_cube.data) # 0th freq

    CS = plt.contourf(W, F, pow_cube.data, levels=levels,
                      cmap=getcolors('WhBlGrYeRe'), extend="both")
    plt.colorbar(CS)

    # set axes range
    plt.ylim(0, NW)
    plt.xlim(fMin, fMax)

    # Line markers of periods
    frqs = [80, 30]
    for frq in frqs:
        plt.plot([1. / frq, 1. / frq], [-0, 15], 'k--', lw=0.5)
        plt.text(1. / frq, 5.5, str(frq) + 'd', {'color' : 'k'})

    plt.plot([0, 0], [-0, 15], 'k:', lw=0.25)
    plt.title(title)  #
    plt.xlabel('Westward     Frequency     Eastward')
    plt.ylabel('Zonal wavenumber')
    save_figure(figname)#,bbox_inches='tight')
    #plt.show()
    plt.close()

def rmm_plot(rmmfile, figname):

    x = 4
    y = 0.707107
    linewidth = 0.25
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal', facecolor='lightgrey')
    plt.plot(np.arange(10, 11))
    plt.plot([-x, -y], [-x, -y], 'k', lw=linewidth)
    plt.plot([y, x], [y, x], 'k', lw=linewidth)
    plt.plot([-x, -y], [x, y], 'k', lw=linewidth)
    plt.plot([y, x], [-y, -x], 'k', lw=linewidth)
    plt.plot([-x, -1], [0, 0], 'k', lw=linewidth)
    plt.plot([1, x], [0, 0], 'k', lw=linewidth)
    plt.plot([0, 0], [-x, -1], 'k', lw=linewidth)
    plt.plot([0, 0], [1, x], 'k', lw=linewidth)

    c = mpatches.Circle((0, 0), 1, fc="lightgrey", ec="k", lw=linewidth)
    ax.add_patch(c)

    plt.text(-3.75, -2, 'Phase 1', va='center', rotation=90, fontsize='smaller')
    plt.text(-2, -3.75, 'Phase 2', ha='center', fontsize='smaller')
    plt.text(2, -3.75, 'Phase 3', ha='center', fontsize='smaller')
    plt.text(3.75, -2, 'Phase 4', ha='center', va='center', rotation=90, fontsize='smaller')
    plt.text(3.75, 2, 'Phase 5', ha='center', va='center', rotation=90, fontsize='smaller')
    plt.text(2, 3.75, 'Phase 6', ha='center', va='center', fontsize='smaller')
    plt.text(-2, 3.75, 'Phase 7', ha='center', va='center', fontsize='smaller')
    plt.text(-3.75, 2, 'Phase 8', va='center', rotation=90, fontsize='smaller')

    plt.text(-4.75, 0, 'Western Hem, Africa', va='center', rotation=90)
    plt.text(0, -4.75, 'Indian Ocean', ha='center', rotation=0)
    plt.text(4.5, 0, 'Maritime continent', ha='center', va='center', rotation=90)
    plt.text(0, 4.5, 'Western Pacific', va='center', ha='center', rotation=0)

    plt.ylim(-x, x)
    plt.xlim(-x, x)

    if os.path.isfile(rmmfile):
        C = np.loadtxt(rmmfile)
        year = C[:, 0]
        pc1 = C[:, 3]
        pc2 = C[:, 4]
        pha = C[:, 6]

        inds = np.where(pc1 != -999.)[0]
        inds = inds[-365:]  # use max 60 days

        plt.plot(pc1[inds], pc2[inds], 'b')



    save_figure(figname)#,bbox_inches='tight')
    #plt.title('MJO Phase space')#
    #plt.xlabel('RMM1')
    #plt.ylabel('RMM2')
    #plt.show()
    plt.close()

def phasefreq_polar_plot(rmmfile, figname, label):
    df = pd.read_csv(rmmfile, sep=' ', names=['year', 'month', 'day', 'rmm1', 'rmm2', 'phase', 'amp'])
    df['date'] = [y * 10000 + m * 100 + d for y, m, d in zip(df.year, df.month, df.day)]
    df = df.loc[df.phase != -999]

    # Compute probailities
    probs = []
    probs.append(len(df.loc[df.amp < 1.]) / len(df))
    for phase in range(1, 9):
        probs.append(len(df.loc[(df.amp >= 1) & (df.phase == phase)]) / len(df))

    probs_df = pd.DataFrame([(i, p * 100) for i, p in enumerate(probs)], columns=['phase', 'prob'])

    _set_plot_theme()
    # Compute pie slices
    # Just pull the 8 phases for plotting
    eight_phases = probs_df.loc[probs_df.phase > 0]
    num_phases = len(eight_phases)
    angle_start = -157.5
    eight_phases['angle'] = np.linspace(angle_start, -angle_start, num_phases) * np.pi / 180
    eight_phases['radii'] = eight_phases["prob"] * 1.
    eight_phases['width'] = 1.5 * np.pi / 4 * eight_phases['radii'] / 30.
    my_cmap = plt.cm.get_cmap('GnBu')
    # radii_normalized = (eight_phases['radii'] - eight_phases['radii'].min()) / (eight_phases['radii'].max() - eight_phases['radii'].min())
    radii_normalized = (eight_phases['radii'] - 6) / 2.
    eight_phases['data_color'] = [my_cmap(value) for value in radii_normalized]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, aspect='equal', projection='polar')  # , axisbg='lightgrey')
    cs = ax.bar(eight_phases['angle'], eight_phases['radii'], width=eight_phases['width'], bottom=1,
                color=eight_phases['data_color'], alpha=0.7, edgecolor='grey')
    ax.set_ylim(-1, 11)
    ax.yaxis.grid(True)
    ax.xaxis.grid(True)
    ax.set_xticklabels([])
    ax.set_rlabel_position(0)  # Set the radial labels at 45 degrees
    ax.patch.set_alpha(1)
    # Customize the plot
    ax.set_theta_direction(1)  # Clockwise direction for the phases

    circle = Circle((0, 0), 1.9, transform=ax.transData._b, color='m', alpha=1, fill=True)
    ax.add_patch(circle)

    for i, row in eight_phases.iterrows():
        print(i, row['prob'], row['angle'])
        ax.text(row['angle'], row['radii'], f"{row['prob']:.2f}%", fontsize=10,
                horizontalalignment='center', verticalalignment='center')
        ax.text(row['angle'], 10.5, f"P{str(row['phase'])}", horizontalalignment='center', verticalalignment='center')
    no_mjo_prob = f"Weak \n No MJO \n{probs_df.loc[probs_df.phase == 0]['prob'].values[0]:.2f}%"
    ax.text(0, -1, no_mjo_prob, horizontalalignment='center', verticalalignment='center', fontsize=12)
    title = f'MJO Frequency \n{label}'
    if len(title) > 40:
        title = textwrap.fill(title, width=40, break_long_words=False,
                                break_on_hyphens=False)
    ax.set_title(title, fontsize=14)
    # ax.plot(eight_phases['angle'], eight_phases['radii'])
    rad = 10
    ax.text(-np.pi, rad, 'Western Hem, Africa', va='center', rotation=90)
    ax.text(-np.pi / 2, rad, 'Indian Ocean', ha='center', rotation=0)
    ax.text(0, rad, 'Maritime continent', ha='center', va='center', rotation=90)
    ax.text(np.pi / 2, rad, 'Western Pacific', va='center', ha='center', rotation=0)

    save_figure(figname)  # ,bbox_inches='tight')
    # plt.show()
    plt.close()
    print('Plotted %s' %figname)

def phasefreq_plot(rmmfile, figname):
    print(rmmfile)
    if os.path.isfile(rmmfile):
        C = np.loadtxt(rmmfile)
        year = C[:, 0]
        month = C[:, 1]
        pc1 = C[:, 3]
        pc2 = C[:, 4]
        amp = C[:, 5]
        pha = C[:, 6]

        freq = np.zeros((2, 9)) # 12 months, 8 phases
        for p in np.arange(8) + 1:
            # Summer
            freq[0, p] = len([i for i in range(len(pha)) if  amp[i] > 1.0 and pha[i] == p and (month[i] >= 5 and month[i] <= 10)])
            freq[0, 0] = len([i for i in range(len(pha)) if  amp[i] <= 1.0 and (month[i] >= 5 and month[i] <= 10) ])
            # Winter
            freq[1, p] = len([i for i in range(len(pha)) if  amp[i] > 1.0 and pha[i] == p and (month[i] >= 11 or month[i] <= 4) ])
            freq[1, 0] = len([i for i in range(len(pha)) if  amp[i] <= 1.0 and (month[i] >= 11 or month[i] <= 4) ])


        freq = 100.*freq / freq.sum()
        phases = np.arange(0, 9, 1)
        months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        fig = plt.figure(1)
        ax = fig.add_subplot(111, aspect=0.9, facecolor='lightgrey')

        p1 = plt.bar(phases[1:] - 0.35, freq[0, 1:], color='w', width=0.35)
        p2 = plt.bar(phases[1:] + 0.0, freq[1, 1:], color='g', width=0.35)

        plt.text (1.5, 4.25, 'No MJO summer : ' + str(round(freq[0, 0])))
        plt.text (1.5, 4., 'No MJO winter : ' + str(round(freq[1, 0])))
        ax.yaxis.grid(color='lightgray', linestyle='-')
        plt.ylim([0, 5])
        plt.legend((p1[0], p2[0]), ('Summer', 'Winter'))
        plt.xlabel('MJO phases')
        plt.ylabel('% freqency')
        plt.title('')
        save_figure(figname)#,bbox_inches='tight')
        #plt.show()
        plt.close()

