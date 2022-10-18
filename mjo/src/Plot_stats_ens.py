#!/usr/local/sci/bin/python2.7
import matplotlib.pyplot as plt
import matplotlib 
import numpy as np
import matplotlib.colors as mcolors
import model_info
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rc('font', family='Times New Roman') 

def get_named_colors(colors, sort_colors=True, emptycols=0):

    # Sort colors by hue, saturation, value and name.
    if sort_colors is True:
        by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(color))),
                         name)
                        for name, color in colors.items())
        names = [name for hsv, name in by_hsv]
    else:
        names = list(colors)
    return names

def plot_cases_hist(runids, plotname='Plot_stats.ps'):
    cases = ['ALL', 'IO', 'MC', 'WP']
    
    plt.figure(figsize=(10,10), dpi=100)
    
    bar_width = 0.15

    #colors = ['royalblue', 'forestgreen', 'darkgoldenrod', 'darkorchid', 'lavender', 'lightgreen']
    # Get Tableau colors
    # https://matplotlib.org/gallery/color/named_colors.html#sphx-glr-gallery-color-named-colors-py
    colors = get_named_colors(mcolors.TABLEAU_COLORS, sort_colors=False, emptycols=0)
    print(colors)
    stats = np.zeros((len(cases), 4)) 
    for r, runid in enumerate(runids):
        filename = "STATS_ENS_%s_MJO.txt" % (runid)
        print(filename)
        f = open(filename, "r")
        lines = f.readlines()
        for c, line in enumerate(lines):
            stats[c, :] = [float(x) for x in line.split()[:]]
        
        print(stats)
        yrange = [2*rx + r * bar_width for rx in range(len(cases))]
        print(yrange)
        xmean = 100.*stats[:,2]/stats[0,2]
        x25   = 100.*stats[:,0]/stats[0,2]
        x75   = 100.*stats[:,3]/stats[0,2]
        print(xmean, x25, x75)
        
        #plt.barh(yrange[1:], xmean[1:], bar_width, color=colors[r])
        plt.plot(xmean[1:], yrange[1:], 'o', color=colors[r], markersize=10)
        for c in range(1, len(cases), 1):
            plt.plot([x25[c], x75[c]], [yrange[c], yrange[c]], linewidth=3, color=colors[r])
            plt.text(x75[c]+1, yrange[c], model_info.label(runid), rotation=0, va='center')
    #plt.ylim([0.5, len(cases)])
    plt.xlim([0., 100])
    plt.grid()
    ax = plt.gca()
    txt_xrange = [2*rx + bar_width * len(cases) for rx in range(len(cases))]
    ax.yaxis.grid(False)
    ax.set_yticks(txt_xrange[1:])
    ax.set_yticklabels(cases[1:])

    plt.ylabel('MJO type', fontsize=18)
    plt.xlabel('Percentage of all events', fontsize=18)
    plt.title('MJO events', fontsize=20)
    
    plt.savefig(plotname)
    plt.show()

if __name__ == '__main__':

    runids = ['obs','u-ab680','u-ab674', 'u-ae871', 'be034', 'bd818', 'be402', 'be408', 'be362']
    pltname = 'Plot_stats_ENS_GAGC_GOML_Runs.ps'
    plot_cases_hist(runids, plotname=pltname)
    
    

