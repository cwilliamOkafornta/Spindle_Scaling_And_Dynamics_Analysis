###########################################################################################################################
'''
This script is used as a library for a plot function PlotCellVolumeLength.
The PlotCellVolumeLength function compute for the following plot parameters:
    - pole-to-pole final distance vs cell length
    - pole-to-chromosome final distance vs cell length
    - chromosome-to-chromosome distance vs cell length
The output comprise of one figure showing all the plots stated above. Also, the R² and slope values are saved as a table 
in a csv file.
'''
###########################################################################################################################

# library packages
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import linregress

###########################################################################################################################

def PlotCellVolumeLength(folder_input, folder_output, file_input, plotfilename, savefilename):
    # read into the folder and assign variable to the desired file
    file_list = []
    for file in file_input:
        filepath = os.path.join(folder_input, file)
        try:
            df = pd.read_csv(filepath, index_col=None, encoding='utf-8')
        except:
            df = pd.read_csv(filepath, index_col=None, encoding='latin-1')
        file_list.append(df)

    dfTable = file_list[0]

    # plot axes parameters 
    x = dfTable['Cell length (µm)']
    y_pp = dfTable['P-to-P (µm)']
    y_pc = dfTable['P-to-C (µm)']
    y_cc = dfTable['C-to-C (µm)']

    # plot errors
    x_err = dfTable['Error propagation']
    y_err_pp = dfTable['P-to-P_std']
    y_err_pc = dfTable['P-to-C_std']
    y_err_cc = dfTable['C-to-C_std']

    # define the plot stle and size
    sns.set_style(style=None)
    plt.figure(figsize=(7,8))

    # define the types of marker shapes (see matplotlib marker letters[https://matplotlib.org/stable/api/markers_api.html])
    marker = ['o', 's', 'D', '^', 'H', 'P', '*', 'p']

    '''Pole-to-pole'''
    ax = sns.regplot(data=dfTable, x=x, y=y_pp, scatter=False, line_kws={'color': '#159C6A'}, ci=None)
    for i, (x_val, y_val, marker_val) in enumerate(zip(x, y_pp, marker)):
        plt.scatter(x_val, y_val, marker=marker_val, color='#F1671C', s=50)
    ax.errorbar(x=x, y=y_pp, xerr=x_err, yerr=y_err_pp, capsize=2, elinewidth=1 , fmt='none', color='#F11F1C')
    reg_line_pp = mpatches.Patch(color='#159C6A', label='Pole-to-pole', linewidth=5)

    '''Pole-to-chromosome'''
    ax1 = sns.regplot(data=dfTable, x=x, y=y_pc, scatter=False, line_kws={'color': '#1C3FF1'}, ci=None, ax=ax)
    for i, (x_val, y_val, marker_val) in enumerate(zip(x, y_pc, marker)):
        plt.scatter(x_val, y_val, marker=marker_val, color='#F1671C', s=50)
    ax1.errorbar(x=x, y=y_pc, xerr=x_err, yerr=y_err_pc, capsize=2, elinewidth=1 , fmt='none', color='#F11F1C')
    reg_line_pc = mpatches.Patch(color='#1C3FF1', label='Pole-to-chromosome', linewidth=5)

    '''Chromosome-to-chromosome'''
    ax2 = sns.regplot(data=dfTable, x=x, y=y_cc, scatter_kws={'s': 20}, line_kws={'color': '#F11CAA'}, ci=None, ax=ax, label='Chromosome-to-chromosome')
    for i, (x_val, y_val, marker_val) in enumerate(zip(x, y_cc, marker)):
        plt.scatter(x_val, y_val, marker=marker_val, color='#F1671C', s=50)
    ax2.errorbar(x=x, y=y_cc, xerr=x_err, yerr=y_err_cc, capsize=2, elinewidth=1 , fmt='none', color='#F11F1C')
    reg_line_cc = mpatches.Patch(color='#F11CAA', label='Chromosome-to-chromosome', linewidth=5)

    # plot properties
    plt.xlim(30, 5)
    plt.ylim(0, 30)
    plt.tick_params(axis='x', labelsize=20)
    plt.tick_params(axis='y', labelsize=20)
    # plt.title('Final distance vs Cell length', fontsize=15)
    plt.xlabel('Cell length [µm]', fontsize=30)
    plt.ylabel('Final distance [µm]', fontsize=30)
    plt.legend(handles=[reg_line_pp, reg_line_pc, reg_line_cc])
    plt.savefig(os.path.join(folder_output, plotfilename+'.png'), dpi=600)
    plt.savefig(os.path.join(folder_output, plotfilename+'.svg'), dpi=600)

    # linear regression coefficients
    '''Pole-to-pole'''
    slope_pp, intercept_pp, r_value_pp, p_value_pp, std_err_pp = linregress(x, y_pp)

    '''Pole-to-chromosome'''
    slope_pc, intercept_pc, r_value_pc, p_value_pc, std_err_pc = linregress(x, y_pc)

    '''Chromosome-to-chromosome'''
    slope_cc, intercept_cc, r_value_cc, p_value_cc, std_err_cc = linregress(x, y_cc)

    # R² values
    table_r_value = {'Pole-to-pole vs cell length': r_value_pp, 'Pole-to-chromosome vs cell length': r_value_pc, 'Chromosome-to-chromosome vs cell length': r_value_cc}
    rTable = pd.DataFrame.from_dict(table_r_value, orient='index', columns=['R² value'])

    # slope
    table_slope = {'Pole-to-pole vs cell length': slope_pp, 'Pole-to-chromosome vs cell length': slope_pc, 'Chromosome-to-chromosome vs cell length': slope_cc}
    slopeTable = pd.DataFrame.from_dict(table_slope, orient='index', columns=['Slope'])

    # merge dataframes
    valueTable = pd.concat([slopeTable, rTable], axis=1)
    valueTable.index.name = 'Parameter'

    # save as csv
    valueTable.to_csv(os.path.join(folder_output, savefilename+'.csv'))
    
    plt.close()