####################################################################################################################################################
'''
This script is used to calculate for the following parameters:
    - cell length from the cell volume
    - cell volume propagation error (to be used for the error bars for the cell length)
    - mean values for each cell stage of the spindle length measurement
    - standard deviation from the calculated mean values.
The calculated parameters are appended in a new table and exported as a csv file to be used for plotting and further analysis.
'''
####################################################################################################################################################

# library packages
import os
import pandas as pd
import numpy as np
####################################################################################################################################################

def cellVolumeLengthValues(folder_input, folder_output, file_input, savefilename):
    # read out the cvs files
    file_list = []
    for file in file_input:
        file_path = os.path.join(folder_input, file)
        try:
            df = pd.read_csv(file_path, index_col=None, encoding='utf-8')
        except:
            df = pd.read_csv(file_path, index_col=None, encoding='latin-1')
        file_list.append(df)

    # assign new variable to each dataFrame
    cell_vol = file_list[0]
    chromosomes_dist = file_list[1]
    poles_dist = file_list[2]
    pole_chromosomes = file_list[3]

    # calculate for cell length by calculating the cubic square of the cell volume
    cell_length = cell_vol['Vol (µm³)']**(1/3)
    cell_vol['Cell length [µm]'] = cell_length

    # calculate for the propagation error of uncertainty of the cell length
    '''
    Propagation error of uncertainty for the cell volume is calculated using the following expression
        propagation error of cell length = |(1/3)Volume^(-2/3) * standard deviation of volume|
    '''
    Error_prop = abs((1/3) * cell_vol['Vol (µm³)']**(-2/3) * cell_vol['std_Vol'])
    cell_vol['Error propagation'] = Error_prop

    # chromosome-to-chromosome final distance
    dist_mean_chromosomes = pd.DataFrame()
    mean_chromosomes_std = pd.DataFrame()
    for col in chromosomes_dist.columns:
        '''mean chromosome-to-chromosome final distance'''
        new_val_temp = pd.DataFrame.from_dict({f"{col}": np.mean(chromosomes_dist[col])}, orient='index', columns=['C-to-C (µm)'])
        dist_mean_chromosomes = pd.concat([dist_mean_chromosomes, new_val_temp], axis=0)

        '''standard deviation chromosome-to-chromosome final distance'''
        new_chromosomes_std = pd.DataFrame.from_dict({f"{col}": np.std(chromosomes_dist[col])}, orient='index', columns=['C-to-C_std']) 
        mean_chromosomes_std = pd.concat([mean_chromosomes_std, new_chromosomes_std], axis=0)

    # pole-to-pole final distance
    dist_mean_poles = pd.DataFrame()
    mean_poles_std = pd.DataFrame()
    for col_poles in poles_dist.columns:
        '''mean pole-to-pole final distance'''
        new_val_pole_temp = pd.DataFrame.from_dict({f"{col_poles}": np.mean(poles_dist[col_poles])}, orient='index', columns=['P-to-P (µm)'])
        dist_mean_poles = pd.concat([dist_mean_poles, new_val_pole_temp], axis=0)

        '''standard deviation pole-to-pole final distance'''
        new_pole_std = pd.DataFrame.from_dict({f"{col_poles}": np.std(poles_dist[col_poles])}, orient='index', columns=['P-to-P_std'])
        mean_poles_std = pd.concat([mean_poles_std, new_pole_std], axis=0)

    # pole-to-chromosome final distance 
    dist_mean_PC = pd.DataFrame()
    mean_PC_std = pd.DataFrame()
    for col_PC in pole_chromosomes:
        '''mean pole-to-chromosome final distance'''
        new_val_PC_temp = pd.DataFrame.from_dict({f"{col_PC}": np.mean(pole_chromosomes[col_PC])}, orient='index', columns=['P-to-C (µm)'])
        dist_mean_PC = pd.concat([dist_mean_PC, new_val_PC_temp], axis=0)

        '''standard deviation pole-to-chromosome final distance'''
        new_val_PC_std = pd.DataFrame.from_dict({f"{col_PC}": np.std(pole_chromosomes[col_PC])}, orient='index', columns=['P-to-C_std'])
        mean_PC_std = pd.concat([mean_PC_std, new_val_PC_std], axis=0)

    newTable = pd.concat([dist_mean_chromosomes, mean_chromosomes_std, dist_mean_poles, mean_poles_std, dist_mean_PC, mean_PC_std], axis=1).reset_index(drop=True)
    newTable = pd.concat([cell_vol, newTable], axis=1)
    newTable.to_csv(os.path.join(folder_output, savefilename+'.csv'), index=False) # save csv