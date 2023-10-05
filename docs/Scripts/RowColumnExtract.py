### Row and column extraction script ###

''''
This is script is used to extract the X, Y, and Z coordinates of the 2 spindle poles.
The input csv file with table is generated after running a spot detection script to track the 
2 poles. The wrongly detected spots are corrected on the ROI Manager in Fiji and a csv file
containing the coordinate information is saved, which will subsequently serve as the input file 
for this script.
'''

# library packages
import pandas as pd
import numpy as np
import os

# define a function for the operation.
def SpotCoordinateExtract(folder, save_folder, input_file, output_file):

    # create a list for the new table
    newTable = []

    # read out the csv file
    df = pd.read_csv(os.path.join(folder, input_file))

    # filter the DataFrame to select only First_Spot and Second_Spot rows
    spot_df = df[df['Name'].str.contains('First_Spot|Second_Spot')]

    # select the X, Y, and Z columns for each First_Spot and Second_Spot row
    first_spot_df = spot_df[spot_df['Name'].str.startswith('First_Spot')][['X', 'Y', 'Z']]
    second_spot_df = spot_df[spot_df['Name'].str.startswith('Second_Spot')][['X', 'Y', 'Z']]

    # rename the columns to X1, Y1, Z1 and X2, Y2, Z2
    first_spot_df.columns = ['X1', 'Y1', 'Z1']
    second_spot_df.columns = ['X2', 'Y2', 'Z2']

    # reset the index column
    first_spot_df = first_spot_df.reset_index(drop=True)
    second_spot_df = second_spot_df.reset_index(drop=True)

    # merge the 2 dfs
    merged_df = pd.concat([first_spot_df, second_spot_df], axis=1)
    merged_df['frame'] = range(len(merged_df))  # add the new frame column

    # saveTable = newTable.append(merged_df)
    newTable += [merged_df]

    merged_table = pd.concat(newTable, ignore_index=True)
    merged_table = merged_table.rename_axis('frame')
    merged_table.to_csv(os.path.join(save_folder, output_file), index=False)
