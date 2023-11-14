# library packages for analyses
import os 
import csv
import math
import glob
import pandas as pd
import statistics
import numpy as np
from numpy import mean

######################################################################################################################################
#                                                                                                                                    #
# The following functions are used to process the csv files generated after  extracting the coordinate positions information of the  #
# centrosomes and chromosomes intensities from kymographs. The coordinate positions information can be extracted with a MatLab or    #
# python custom code (see folder: intensity_readout_scripts).                                                                        #
# This script is tested for the csv files generated from MatLab code.                                                                #
#                                                                                                                                    #
# The workflow is adopted to give a new heading to each column on the csv files, followed by calculating the distances using the     #
# positions provided in each file. Each file is labeled with the name of the kymograph from a given cell type of the C. elegans      #
# embryo.                                                                                                                            #
#                                                                                                                                    #
# The next step is to extract similar columns in each file and add them to a created csv file with a given name. The new files are   #
# saved in different folder.                                                                                                         #
#                                                                                                                                    #
######################################################################################################################################

def distanceCalculation(files):
    for file in files:
        try:
            # define the title of each header of each columns with values
            title_head = ['frame', 'pole1', 'pole2', 'chromatid1', 'chromatid2', 'No. of chromosome']

            # read out the csv file 
            df = pd.read_csv(file, index_col=None, names=title_head)

            # calculate the time in seconds from the frame column
            df['time'] = (df['frame']-11)*10.3

            '''calculate the distances by parsing the values from the given dataframe (columns)'''
            list_pole_chromosome = [] # create a list of the pole-to-chromosome distances of both sides

            # calculate the 3D Euclidean distance of poles and chromosomes
            p1_p2 = ((abs(df['pole1']-df['pole2']))*0.1)**2
            p1_c1 = ((abs(df['pole1']-df['chromatid1']))*0.1)**2
            p2_c2 = ((abs(df['pole2']-df['chromatid2']))*0.1)**2
            c1_c2 = ((abs(df['chromatid1']-df['chromatid2']))*0.1)**2

            df['pole1_pole2'] = round(((p1_p2)**0.5), 4)
            df['pole1_chromatid1'] = round(((p1_c1)**0.5), 4)
            df['pole2_chromatid2'] = round(((p2_c2)**0.5), 4)
            df['chromatid1_chromatid2'] = round(((c1_c2)**0.5), 4)

            # calculate the mean of the pole_to_chromosome distances (for the 2 poles)
            list_pole_chromosome.append(df['pole1_chromatid1'])
            list_pole_chromosome.append(df['pole2_chromatid2'])
            chromosome_pole = list_pole_chromosome
            mean_pole_chromosome = [mean(i) for i in zip(*chromosome_pole)]

            # standard devaition using the P [population] factor as a test input, i.e (N-1) 
            std_pole_chromosome = [np.std(i) for i in zip(*chromosome_pole)]

            df['Avg_pole_chromosome'] = mean_pole_chromosome
            df['std_pole_chromosome'] = std_pole_chromosome

            # save files in the same folder with the same filename
            df.to_csv(os.path.join(r'', file), index=False)
        
        except Exception:
            pass
        
    return None

#####################################################################################################################################
# Both "parameterExtract()" and "extractColumn()" have the same function                                                            #
#####################################################################################################################################

def parameterExtract(file_path, folder_input, folder_output):
    parameters = ['pole1_pole2', 'pole1_chromatid1', 'pole2_chromatid2', 'chromatid1_chromatid2',
                  'Avg_pole_chromosome', 'std_pole_chromosome']
    
    ''' create an empty dictionary and then extract the columns with a given header and add to the dictionary'''
    # create a dictionary to add the list of dataframes 
    newFiles = {}
    for parameter in parameters:
        newFiles[parameter] = pd.DataFrame() # create an empty dataframe for the parameter list in the dictionary
        
        # read out the files in the folder
        with os.scandir(file_path) as folder_input:
            for j in folder_input:
                df_files = pd.read_csv(j, index_col=None)

                # assign the extracted columns with the parent filename
                try:
                    filename = j.name
                    columnname = filename.split('.')[0]
                    newFiles[parameter][columnname] = df_files[parameter]
                except Exception:
                    pass
        
        # save the individual table with similar extracted columns in a different folder
        newFiles[parameter].to_csv(os.path.join(folder_output, parameter+'.csv'), index=False)

#####################################################################################################################################
        
def extractColumn(folder_input, folder_output):
    '''
    Create a list of all the new files
    '''
    Exp_Avg_pole_chromosome = []
    Exp_chromatid_chromatid = []
    Exp_pole_pole = []
    Exp_std_pole_chromosome = []

    with os.scandir(folder_input) as files_1:
        # read out the csv file 
        for i in files_1: 
            df= pd.read_csv(i, index_col=None)

            # Average pole-to-chromosome distance
            value_column_p_c = {i: (df['Avg_pole_chromosome'])}
            value_measure_p_c = pd.DataFrame(value_column_p_c)
            value_rename_p_c = value_measure_p_c.rename(columns={i:i.name.split('.')[0]}) # remove the file ext/os.scandir attribute
            Exp_Avg_pole_chromosome.append(value_rename_p_c)

            # chromosome-to-chromosome distance
            value_column_c_c = {i: (df['chromatid1_chromatid2'])}
            value_measure_c_c = pd.DataFrame(value_column_c_c)
            value_rename_c_c = value_measure_c_c.rename(columns={i:i.name.split('.')[0]}) # remove the file ext/os.scandir attribute
            Exp_chromatid_chromatid.append(value_rename_c_c)

            # pole-to-pole distance
            value_column_p_p = {i: (df['pole1_pole2'])}
            value_measure_p_p = pd.DataFrame(value_column_p_p)
            value_rename_p_p = value_measure_p_p.rename(columns={i:i.name.split('.')[0]}) # remove the file ext/os.scandir attribute
            Exp_pole_pole.append(value_rename_p_p)

            # standard deviation of pole-to-chromosome distance using the p [population] factor, i.e (N-1)
            value_column_std_p_c = {i: (df['std_pole_chromosome'])}
            value_measure_std_p_c = pd.DataFrame(value_column_std_p_c)
            value_rename_std_p_c = value_measure_std_p_c.rename(columns={i:i.name.split('.')[0]}) #remove the file ext/os.scandir attribute
            Exp_std_pole_chromosome.append(value_rename_std_p_c)
    
    # combine similar columns into a new table
    Exp_combine_p_c = pd.concat(Exp_Avg_pole_chromosome, ignore_index=False, axis=1)
    Exp_combine_c_c = pd.concat(Exp_chromatid_chromatid, ignore_index=False, axis=1)
    Exp_combine_p_p = pd.concat(Exp_pole_pole, ignore_index=False, axis=1)
    Exp_combine_std_p_c = pd.concat(Exp_std_pole_chromosome, ignore_index=False, axis=1)
    
    # save the new created tables with given filenames
    p_c = Exp_combine_p_c.to_csv(os.path.join(folder_output, 'Avg_pole_chromosome.csv'),index=False)
    c_c = Exp_combine_c_c.to_csv(os.path.join(folder_output, 'chromatid_chromatid.csv'), index=False)
    p_p = Exp_combine_p_p.to_csv(os.path.join(folder_output, 'pole_pole.csv'), index=False)
    std_p_c = Exp_combine_std_p_c.to_csv(os.path.join(folder_output, 'std_pole_chromosome.csv'), index=False)

    return p_c, c_c, p_p, std_p_c

#####################################################################################################################################

