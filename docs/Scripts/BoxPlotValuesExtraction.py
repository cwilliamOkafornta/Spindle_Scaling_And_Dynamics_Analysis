import os
import pandas as pd

###############################################################################################################################
'''
The functions below are used for extracting multiple column values from a csv file and creating a new table where the values
of different cell types are grouped by their cell stages. The goal is to create a new DataFrame that will allow for making
an understandable boxplot for the length measurement based on the cell stages.
'''

###############################################################################################################################

''' For Pole_to_pole distance extraction'''

# dictionary for each table
new_table_initial = {}
new_table_final = {}
new_table_elongation = {}
new_table_Metaphase_len = {}
new_table_reduction = {}
new_table_ratio = {}
new_table_speed = {}

def PoleToPole_CSVvalues_BoxPlot(input_folder, output_folder):

    # read out files from the folder
    with os.scandir(input_folder) as file_path:
        for file in file_path:
            if file.name.endswith('.csv'):
                Result_fit_table = pd.read_csv(file.path, index_col=None)

                # Extract the cell types from the 'Cells' column
                Result_fit_table['Cell Type'] = Result_fit_table['Cells'].str.split('_', expand=True)[0]

                # iterate over the Cell Type column
                for cell_type in Result_fit_table['Cell Type'].unique():

                    # Initial length
                    cell_type_column_initial = f'{cell_type}'
                    values_initial = Result_fit_table .loc[Result_fit_table ['Cell Type'] == cell_type, 
                                                           'Initial pole_pole length (µm)'].reset_index(drop=True)
                    new_table_initial[cell_type_column_initial] = values_initial

                    # final length
                    cell_type_column_final = f'{cell_type}'
                    values_final = Result_fit_table .loc[Result_fit_table['Cell Type'] == cell_type, 
                                                         'Final pole_pole length (µm)'].reset_index(drop=True)
                    new_table_final[cell_type_column_final] = values_final

                    # elongation rate
                    cell_type_column_elongation = f'{cell_type}'
                    values_elongation = Result_fit_table .loc[Result_fit_table['Cell Type'] == cell_type, 
                                                              'Elongation rate (µm/min)'].reset_index(drop=True)
                    new_table_elongation[cell_type_column_elongation] = values_elongation

                    # metaphase pole-to-pole length
                    cell_type_column_metaphase = f'{cell_type}'
                    values_metaphase = Result_fit_table .loc[Result_fit_table['Cell Type'] == cell_type, 
                                                              'Metaphase_length (µm/min)'].reset_index(drop=True)
                    new_table_Metaphase_len[cell_type_column_metaphase] = values_metaphase
                    
                    # Create the new table DataFrame
                    newTable_initial = pd.DataFrame(new_table_initial)
                    newTable_final = pd.DataFrame(new_table_final)
                    newTable_elongation = pd.DataFrame(new_table_elongation)
                    newTable_metaphase = pd.DataFrame(new_table_Metaphase_len)
                    
    # save in another folder
    newTable_initial.to_csv(os.path.join(output_folder, 'Initial pole_to_pole_length ().csv'), index=False, encoding='utf-8')
    newTable_final.to_csv(os.path.join(output_folder, 'Final pole_pole length ().csv'), index=False, encoding='utf-8')
    newTable_elongation.to_csv(os.path.join(output_folder, 'Elongation rate pole_to_pole ().csv'), index=False, encoding='utf-8')
    newTable_metaphase.to_csv(os.path.join(output_folder, 'Metaphase_pole_to_pole_length ().csv'), index=False, encoding='utf-8')
    
###############################################################################################################################

''' For Pole_to_Chromosome distance extraction'''

def PoleToChromosome_CSVvalues_BoxPlot(input_folder, output_folder):
    # read out files from the folder
    with os.scandir(input_folder) as file_path:
        for file in file_path:
            if file.name.endswith('.csv'):
                Result_fit_table = pd.read_csv(file.path, index_col=None)

                # Extract the cell types from the 'Cells' column
                Result_fit_table['Cell Type'] = Result_fit_table['Cells'].str.split('_', expand=True)[0]

                # iterate over the Cell Type column
                for cell_type in Result_fit_table['Cell Type'].unique():
                    
                    # Initial length
                    cell_type_column_initial = f'{cell_type}'
                    values_initial = Result_fit_table .loc[Result_fit_table ['Cell Type'] == cell_type, 
                                                           'Initial length (µm)'].reset_index(drop=True)
                    new_table_initial[cell_type_column_initial] = values_initial
                    
                    # final length
                    cell_type_column_final = f'{cell_type}'
                    values_final = Result_fit_table .loc[Result_fit_table['Cell Type'] == cell_type, 
                                                         'Final length (µm)'].reset_index(drop=True)
                    new_table_final[cell_type_column_final] = values_final
                    
                    # rate of distance reduction
                    cell_type_column_reduction = f'{cell_type}'
                    values_reduction = Result_fit_table .loc[Result_fit_table['Cell Type'] == cell_type, 
                                                              'Rate of distance reduction (µm/min)'].reset_index(drop=True)
                    new_table_reduction[cell_type_column_reduction] = values_reduction
                    
                    # ratio (final length vs initial length)
                    cell_type_column_ratio = f'{cell_type}'
                    values_ratio = Result_fit_table .loc[Result_fit_table['Cell Type'] == cell_type, 
                                                              'Ratio (final length / initial length)'].reset_index(drop=True)
                    new_table_ratio[cell_type_column_ratio] = values_ratio
                    
                    # metaphase pole-to-pole length
                    cell_type_column_metaphase = f'{cell_type}'
                    values_metaphase = Result_fit_table .loc[Result_fit_table['Cell Type'] == cell_type, 
                                                              'Metaphase_length (µm/min)'].reset_index(drop=True)
                    new_table_Metaphase_len[cell_type_column_metaphase] = values_metaphase
                    
                    # Create the new table DataFrame
                    newTable_initial = pd.DataFrame(new_table_initial)
                    newTable_final = pd.DataFrame(new_table_final)
                    newTable_reduction = pd.DataFrame(new_table_reduction)
                    newTable_ratio = pd.DataFrame(new_table_ratio)
                    newTable_metaphase = pd.DataFrame(new_table_Metaphase_len)
            
    # save in another folder
    newTable_initial.to_csv(os.path.join(output_folder, 'Initial pole_to_chromosome_length ().csv'), index=False, encoding='utf-8')
    newTable_final.to_csv(os.path.join(output_folder, 'Final pole_to_chromosome_length ().csv'), index=False, encoding='utf-8')
    newTable_reduction.to_csv(os.path.join(output_folder, 'Reduction_rate ().csv'), index=False, encoding='utf-8')
    newTable_ratio.to_csv(os.path.join(output_folder, 'Ratio_final_vs_initial ().csv'), index=False, encoding='utf-8')
    newTable_metaphase.to_csv(os.path.join(output_folder, 'Metaphase_pole_to_chromosome_length ().csv'), index=False, encoding='utf-8')
    

###############################################################################################################################
                    
''' For Chromosome_to_Chromosome distance extraction'''
    
def ChromosomeToChromosome_CSVvalues_BoxPlot(input_folder, output_folder):
    # read out files from the folder
    with os.scandir(input_folder) as file_path:
        for file in file_path:
            if file.name.endswith('.csv'):
                Result_fit_table = pd.read_csv(file.path, index_col=None)

                # Extract the cell types from the 'Cells' column
                Result_fit_table['Cell Type'] = Result_fit_table['Cells'].str.split('_', expand=True)[0]
                
                # iterate over the Cell Type column
                for cell_type in Result_fit_table['Cell Type'].unique():
                    
                    # final length
                    cell_type_column_final = f'{cell_type}'
                    values_final = Result_fit_table .loc[Result_fit_table['Cell Type'] == cell_type, 
                                                         'Final chromosomes length (µm)'].reset_index(drop=True)
                    new_table_final[cell_type_column_final] = values_final
                    
                    # Segregation speed
                    cell_type_column_speed = f'{cell_type}'
                    values_speed = Result_fit_table .loc[Result_fit_table['Cell Type'] == cell_type, 
                                                         'Segregation speed (µm/min)'].reset_index(drop=True)
                    new_table_speed[cell_type_column_speed] = values_speed
            
                    # Create the new table DataFrame
                    newTable_final = pd.DataFrame(new_table_final) 
                    newTable_speed = pd.DataFrame(new_table_speed)
                    
    # save in another folder
    newTable_final.to_csv(os.path.join(output_folder, 'Final chromosome_to_chromosome_length ().csv'), index=False, encoding='utf-8')
    newTable_speed.to_csv(os.path.join(output_folder, 'Segregation speed_C_C ().csv'), index=False, encoding='utf-8')