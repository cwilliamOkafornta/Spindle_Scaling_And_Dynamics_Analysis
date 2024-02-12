# library packages
import os 
import pandas as pd
import numpy as np

######################################################################################################################################
'''
This library script is used to extract values from csv files of different cell types and experimental ID. The extracted values are combined 
based on their embryonic cell stage annotation to form a new DataFrame. Each DataFrame is subsequently saved under the embryonic cell stage ID.
'''
######################################################################################################################################

def extract2cellstage_wt(input_folder, output_folder, file_input, parameter):
    
    new_data = []
    for cell_name in parameter:
        data ={}
        data[cell_name] = pd.DataFrame()
        
        with os.scandir(input_folder) as file_input:
            for i in file_input:
                if i.is_file() and i.name.startswith('2-cell') and i.name.endswith('csv'):
                    df = pd.read_csv(i.path, index_col=None)
                    
                    try:
                        columnname = i.name.split('.')[0]
                        data[cell_name][columnname] = df[cell_name]
                    except Exception:
                        pass
                    
        # new df
        new_table = pd.DataFrame(data[cell_name], index=None)
        new_data.append(new_table)
        
    # concatenate
    df_new_table = pd.concat(new_data, axis=1) 
    
    # add new columns
    df_mean = df_new_table.apply(lambda row: np.mean(row.dropna()), axis=1)
    df_median = df_new_table.apply(lambda row: np.median(row.dropna()), axis=1)
    df_std = df_new_table.apply(lambda row: np.std(row.dropna()), axis=1)
    df_n_value = df_new_table.apply(lambda row: row.count(), axis=1)
    
       # Concatenate all new columns with df_new_table
    new_columns = pd.concat([df_mean.rename('mean'), df_median.rename('median'), df_std.rename('std'), df_n_value.rename('n_value')], axis=1)

    # Concatenate df_new_table with new_columns and 'time' column
    df_new_table = pd.concat([df_new_table, new_columns, df['time']], axis=1)
    
    # save
    df_new_table.to_csv(os.path.join(output_folder, 'WT_2-cell.csv'), index=False)
    
######################################################################################################################################

def extract4cellstage_wt(input_folder, output_folder, file_input, parameter):
    
    new_data = []
    for cell_name in parameter:
        data ={}
        data[cell_name] = pd.DataFrame()
        
        with os.scandir(input_folder) as file_input:
            for i in file_input:
                if i.is_file() and i.name.startswith('4-cell') and i.name.endswith('csv'):
                    df = pd.read_csv(i.path, index_col=None)
                    
                    try:
                        columnname = i.name.split('.')[0]
                        data[cell_name][columnname] = df[cell_name]
                    except Exception:
                        pass
                    
        # new df
        new_table = pd.DataFrame(data[cell_name], index=None)
        new_data.append(new_table)
        
    # concatenate
    df_new_table = pd.concat(new_data, axis=1) 
    
    # add new columns
    df_mean = df_new_table.apply(lambda row: np.mean(row.dropna()), axis=1)
    df_median = df_new_table.apply(lambda row: np.median(row.dropna()), axis=1)
    df_std = df_new_table.apply(lambda row: np.std(row.dropna()), axis=1)
    df_n_value = df_new_table.apply(lambda row: row.count(), axis=1)
    
    # Concatenate all new columns with df_new_table
    new_columns = pd.concat([df_mean.rename('mean'), df_median.rename('median'), df_std.rename('std'), df_n_value.rename('n_value')], axis=1)

    # Concatenate df_new_table with new_columns and 'time' column
    df_new_table = pd.concat([df_new_table, new_columns, df['time']], axis=1)
        
    # save
    df_new_table.to_csv(os.path.join(output_folder, 'WT_4-cell.csv'), index=False)
    
######################################################################################################################################

def extract8cellstage_wt(input_folder, output_folder, file_input, parameter):
    
    new_data = []
    for cell_name in parameter:
        data ={}
        data[cell_name] = pd.DataFrame()
        
        with os.scandir(input_folder) as file_input:
            for i in file_input:
                if i.is_file() and i.name.startswith('8-cell') and i.name.endswith('csv'):
                    df = pd.read_csv(i.path, index_col=None)
                    
                    try:
                        columnname = i.name.split('.')[0]
                        data[cell_name][columnname] = df[cell_name]
                    except Exception:
                        pass
                    
        # new df
        new_table = pd.DataFrame(data[cell_name], index=None)
        new_data.append(new_table)
        
    # concatenate
    df_new_table = pd.concat(new_data, axis=1) 
    
    # add new columns
    df_mean = df_new_table.apply(lambda row: np.mean(row.dropna()), axis=1)
    df_median = df_new_table.apply(lambda row: np.median(row.dropna()), axis=1)
    df_std = df_new_table.apply(lambda row: np.std(row.dropna()), axis=1)
    df_n_value = df_new_table.apply(lambda row: row.count(), axis=1)
    
    # Concatenate all new columns with df_new_table
    new_columns = pd.concat([df_mean.rename('mean'), df_median.rename('median'), df_std.rename('std'), df_n_value.rename('n_value')], axis=1)

    # Concatenate df_new_table with new_columns and 'time' column
    df_new_table = pd.concat([df_new_table, new_columns, df['time']], axis=1)
        
    # save
    df_new_table.to_csv(os.path.join(output_folder, 'WT_8-cell.csv'), index=False)
    
######################################################################################################################################

def extract16cellstage_wt(input_folder, output_folder, file_input, parameter):
    
    new_data = []
    for cell_name in parameter:
        data ={}
        data[cell_name] = pd.DataFrame()
        
        with os.scandir(input_folder) as file_input:
            for i in file_input:
                if i.is_file() and i.name.startswith('16-cell') and i.name.endswith('csv'):
                    df = pd.read_csv(i.path, index_col=None)
                    
                    try:
                        columnname = i.name.split('.')[0]
                        data[cell_name][columnname] = df[cell_name]
                    except Exception:
                        pass
                    
        # new df
        new_table = pd.DataFrame(data[cell_name], index=None)
        new_data.append(new_table)
        
    # concatenate
    df_new_table = pd.concat(new_data, axis=1) 
    
    # add new columns
    df_mean = df_new_table.apply(lambda row: np.mean(row.dropna()), axis=1)
    df_median = df_new_table.apply(lambda row: np.median(row.dropna()), axis=1)
    df_std = df_new_table.apply(lambda row: np.std(row.dropna()), axis=1)
    df_n_value = df_new_table.apply(lambda row: row.count(), axis=1)
    
    # Concatenate all new columns with df_new_table
    new_columns = pd.concat([df_mean.rename('mean'), df_median.rename('median'), df_std.rename('std'), df_n_value.rename('n_value')], axis=1)

    # Concatenate df_new_table with new_columns and 'time' column
    df_new_table = pd.concat([df_new_table, new_columns, df['time']], axis=1)
        
    # save
    df_new_table.to_csv(os.path.join(output_folder, 'WT_16-cell.csv'), index=False)
    
######################################################################################################################################

def extract32cellstage_wt(input_folder, output_folder, file_input, parameter):
    
    new_data = []
    for cell_name in parameter:
        data ={}
        data[cell_name] = pd.DataFrame()
        
        with os.scandir(input_folder) as file_input:
            for i in file_input:
                if i.is_file() and i.name.startswith('32-cell') and i.name.endswith('csv'):
                    df = pd.read_csv(i.path, index_col=None)
                    
                    try:
                        columnname = i.name.split('.')[0]
                        data[cell_name][columnname] = df[cell_name]
                    except Exception:
                        pass
                    
        # new df
        new_table = pd.DataFrame(data[cell_name], index=None)
        new_data.append(new_table)
        
    # concatenate
    df_new_table = pd.concat(new_data, axis=1) 
    
    # add new columns
    df_mean = df_new_table.apply(lambda row: np.mean(row.dropna()), axis=1)
    df_median = df_new_table.apply(lambda row: np.median(row.dropna()), axis=1)
    df_std = df_new_table.apply(lambda row: np.std(row.dropna()), axis=1)
    df_n_value = df_new_table.apply(lambda row: row.count(), axis=1)
    
    # Concatenate all new columns with df_new_table
    new_columns = pd.concat([df_mean.rename('mean'), df_median.rename('median'), df_std.rename('std'), df_n_value.rename('n_value')], axis=1)

    # Concatenate df_new_table with new_columns and 'time' column
    df_new_table = pd.concat([df_new_table, new_columns, df['time']], axis=1)
        
    # save
    df_new_table.to_csv(os.path.join(output_folder, 'WT_32-cell.csv'), index=False)
    
######################################################################################################################################

def extract64cellstage_wt(input_folder, output_folder, file_input, parameter):
    
    new_data = []
    for cell_name in parameter:
        data ={}
        data[cell_name] = pd.DataFrame()
        
        with os.scandir(input_folder) as file_input:
            for i in file_input:
                if i.is_file() and i.name.startswith('64-cell') and i.name.endswith('csv'):
                    df = pd.read_csv(i.path, index_col=None)
                    
                    try:
                        columnname = i.name.split('.')[0]
                        data[cell_name][columnname] = df[cell_name]
                    except Exception:
                        pass
                    
        # new df
        new_table = pd.DataFrame(data[cell_name], index=None)
        new_data.append(new_table)
        
    # concatenate
    df_new_table = pd.concat(new_data, axis=1) 
    
    # add new columns
    df_mean = df_new_table.apply(lambda row: np.mean(row.dropna()), axis=1)
    df_median = df_new_table.apply(lambda row: np.median(row.dropna()), axis=1)
    df_std = df_new_table.apply(lambda row: np.std(row.dropna()), axis=1)
    df_n_value = df_new_table.apply(lambda row: row.count(), axis=1)

    # Concatenate all new columns with df_new_table
    new_columns = pd.concat([df_mean.rename('mean'), df_median.rename('median'), df_std.rename('std'), df_n_value.rename('n_value')], axis=1)

    # Concatenate df_new_table with new_columns and 'time' column
    df_new_table = pd.concat([df_new_table, new_columns, df['time']], axis=1)
        
    # save
    df_new_table.to_csv(os.path.join(output_folder, 'WT_64-cell.csv'), index=False)