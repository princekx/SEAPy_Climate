import os, glob
import pandas as pd
from pandas.api.types import is_numeric_dtype
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from mjo import data_paths
import warnings
warnings.filterwarnings('ignore')
sns.set_palette('colorblind')


def read_model_metrics(area='mjo', aa_home='/home/h03/hadco/public_html/assess/', csv_file='metrics.csv'):
    # read valnote metrics
    csv_file_paths = glob.glob(aa_home + '*/%s/*/%s' % (area, csv_file))
    df_list = []
    for i, csv_file_path in enumerate(csv_file_paths):
        runid = csv_file_path.split('/')[-2]
        dfx = pd.read_csv(csv_file_path).T
        # drop the index row
        dfx = dfx.rename(columns=dfx.iloc[0]).drop(dfx.index[0])
        dfx['runid'] = runid
        df_list.append(dfx)
    df = pd.concat(df_list, axis=0)
    df = df.set_index('runid')
    df = df.drop_duplicates()
    # remove spaces from column names
    df.columns = df.columns.str.replace(' ', '')

    return df


def read_one_model(run, area='mjo', csv_file='metrics.csv'):
    runid = run['runid']
    label = run['label']
    aa_home = os.path.join(data_paths.dirs('data_out_dir'))
    print(aa_home)
    csv_file_paths = glob.glob(aa_home + '/%s/%s' % (runid, csv_file))
    csv_file_path = csv_file_paths[0]
    runid = csv_file_path.split('/')[-2]
    df = pd.read_csv(csv_file_path).T
    df = df.rename(columns=df.iloc[0]).drop(df.index[0])
    df['runid'] = runid
    df = df.set_index('runid')
    # remove spaces from column names
    df.columns = df.columns.str.replace(' ', '')

    return df



def mjo_history(control=None, expt=None):


    out_plot_dir = os.path.join(data_paths.dirs('data_out_dir'), expt['runid'])
    # read all previous model metrics
    df = read_model_metrics()

    # select only obs values
    obs = df.loc['OBS']

    # select only the models
    models = df[~df.index.isin(['OBS'])]
    df_cntl = read_one_model(control)
    df_expt = read_one_model(expt)

    var_names = df_cntl.columns
    # Plot data
    for col_name in var_names:
        print(obs[col_name])
        fig, ax = plt.subplots()
        sns.distplot(models[col_name], ax=ax, label='models')
        sns.distplot(obs[col_name], ax=ax, label='obs')
        yrange = ax.get_ylim()

        yrange_cn = [yrange[0], yrange[1] / 2]
        yrange_ex = [yrange[1] / 2, yrange[1]]
        xrange_cn = [df_cntl[col_name], df_cntl[col_name]]
        xrange_ex = [df_expt[col_name], df_expt[col_name]]

        plt.plot(xrange_cn, yrange_cn, label=control['label'], linewidth=3, alpha=0.7)
        plt.plot(xrange_ex, yrange_ex, label=expt['label'],  linewidth=3, alpha=0.7)
        plt.legend()
        plt.title('')
        # for col_names with '/' in it
        var_name = ''.join(col_name.split('/'))
        figname = os.path.join(out_plot_dir, 'history_%s_ex_%s_vs_con_%s.pdf' % (var_name, expt['runid'], control['runid']))
        plt.savefig(figname, bbox_inches='tight', pad_inches=0)
        plt.close()
        print('Plotted %s' %figname)
