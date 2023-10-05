import os
import numpy as np 
import pandas as pd
from scipy import stats

#######################################################################################################################################################
'''
This script is used to compute for a statistical t-test of unequal variance.
The t-test is computer for the data from different analysis between different conditions of data source; in this case, the data from the wild type and 
RNAi treated conditions of a C. elegans embryo.
The analysis data that this t-test is calculated for are:
    - pole-to-pole distance analysis 
    - pole-to-chromosome distance analysis
    - chromosome-to-chromosome distance analysis
Each function represent the t-test for each distance analysis. 

In each function, the independent t-test is computed for data from similar embryonic stage between two embryonic conditions, and the data between 
two embryonic stages of the same condition.
'''
#######################################################################################################################################################

def ttestPoletoPole(folder_input, filedata_input, folder_output, name_cond1_fileinput, name_cond2_fileinput, condition_1, condition_2):
    # list for files with different conditions
    cond1_final_pp = []
    cond2_final_pp = []
    
    for file in filedata_input:
        filepath = os.path.join(folder_input, file)
        filename = os.path.basename(filepath) # return the name of each file in the folder
        
        # read out the file of condition 1
        if filename == name_cond1_fileinput:
            
            # check for the appropriate format of the csv files
            try:
                df = pd.read_csv(filepath, index_col=None, encoding='utf-8')
            except:
                df = pd.read_csv(filepath, index_col=None, encoding='latin-1')
            cond1_final_pp.append(df)
            
        # read out the file of condition 2
        elif filename == name_cond2_fileinput:
            
            # check for the appropriate format of the csv files 
            try:
                df1 = pd.read_csv(filepath, index_col=None, encoding='utf-8')
            except:
                df1 = pd.read_csv(filepath, index_col=None, encoding='latin-1')
            cond2_final_pp.append(df1)

    # create new dataframes
    cond1_final_pp = pd.concat(cond1_final_pp, ignore_index=True)
    cond2_final_pp = pd.concat(cond2_final_pp, ignore_index=True)
    new_cond1_final_pp = cond1_final_pp.iloc[:,1:] # omit the first column '1-cell'
    
    # compute for Independent t-test between the 2 given conditions
    '''The independent t-test uses statistical test of unequal variance'''
    test_values = pd.DataFrame()
    for column in new_cond1_final_pp.columns:
        cond1_final_pp_column = new_cond1_final_pp[column].dropna()
        cond2_final_pp_column = cond2_final_pp[column].dropna()
        t_stat, t_p_value = stats.ttest_ind(cond1_final_pp_column, cond2_final_pp_column, equal_var=False, nan_policy='omit')
        new_value = {f"{column} of {condition_1} and {condition_2}": t_p_value}
        new_value_df = pd.DataFrame.from_dict(new_value, orient='index', columns=['p-value'])
        test_values = pd.concat([test_values, new_value_df], axis=0)
    test_values.index.name = 'Significance comparison'
    
    # compute for Independet t-test between columns of condition_1 data
    '''The independent t-test uses statistical test of unequal variance'''
    test_own_values_cond1 = pd.DataFrame()
    for i, col1 in enumerate(cond1_final_pp):
        for col2 in cond1_final_pp[i+1:]:
            cond1_final_pp_col1 = cond1_final_pp[col1].dropna()
            cond1_final_pp_col2 = cond1_final_pp[col2].dropna()
            t_stat_1, t_p_value_1 = stats.ttest_ind(cond1_final_pp_col1, cond1_final_pp_col2, equal_var=False, nan_policy='omit')
            own_new_values = {f"{col1} and {col2} of {condition_1}": t_p_value_1}
            own_new_values_df = pd.DataFrame.from_dict(own_new_values, orient='index', columns=['p-value'])
            test_own_values_cond1 = pd.concat([test_own_values_cond1, own_new_values_df], axis=0)
    test_own_values_cond1.index.name = 'Significance comparison'    
    
    # compute for Independet t-test between columns of condition_2 data
    '''The independent t-test uses statistical test of unequal variance'''
    test_own_values_cond2 = pd.DataFrame()
    for i, col1 in enumerate(cond2_final_pp):
        for col2 in cond2_final_pp[i+1:]:
            cond2_final_pp_col1 = cond2_final_pp[col1].dropna()
            cond2_final_pp_col2 = cond2_final_pp[col2].dropna()
            t_stat_2, t_p_value_2 = stats.ttest_ind(cond2_final_pp_col1, cond2_final_pp_col2, equal_var=False, nan_policy='omit')
            own_new_values_cond2 = {f"{col1} and {col2} of {condition_2}": t_p_value_2}
            own_new_values_cond2_df = pd.DataFrame.from_dict(own_new_values_cond2, orient='index', columns=['p-value'])
            test_own_values_cond2 = pd.concat([test_own_values_cond2, own_new_values_cond2_df], axis=0)
    test_own_values_cond2.index.name = 'Significance comparison'
    
    # save dataframes to csv
    test_values.to_csv(os.path.join(folder_output, (f"t-test between {condition_1} and {condition_2} (P_P).csv")))
    test_own_values_cond1.to_csv(os.path.join(folder_output, (f"t-test between embryonic stage of {condition_1} (P_P).csv")))
    test_own_values_cond2.to_csv(os.path.join(folder_output, (f"t-test between embryonic stage of {condition_2} (P_P).csv")))
    
    return None

#######################################################################################################################################################

def ttestChromosometoChromosome(folder_input, filedata_input, folder_output, name_cond1_fileinput, name_cond2_fileinput, condition_1, condition_2):
    # list for files with different conditions
    cond1_final_cc = []
    cond2_final_cc = []
    
    for file in filedata_input:
        filepath = os.path.join(folder_input, file)
        filename = os.path.basename(filepath) # return the name of each file in the folder
        
        # read out the file of condition 1
        if filename == name_cond1_fileinput:
            
            # check for the appropriate format of the csv files
            try:
                df = pd.read_csv(filepath, index_col=None, encoding='utf-8')
            except:
                df = pd.read_csv(filepath, index_col=None, encoding='latin-1')
            cond1_final_cc.append(df)
            
        # read out the file of condition 2
        elif filename == name_cond2_fileinput:
            
            # check for the appropriate format of the csv files 
            try:
                df1 = pd.read_csv(filepath, index_col=None, encoding='utf-8')
            except:
                df1 = pd.read_csv(filepath, index_col=None, encoding='latin-1')
            cond2_final_cc.append(df1)

    # create new dataframes
    cond1_final_cc = pd.concat(cond1_final_cc, ignore_index=True)
    cond2_final_cc = pd.concat(cond2_final_cc, ignore_index=True)
    new_cond1_final_cc = cond1_final_cc.iloc[:,1:] # omit the first column '1-cell'
    
    # compute for Independent t-test between the 2 given conditions
    '''The independent t-test uses statistical test of unequal variance'''
    test_values = pd.DataFrame()
    for column in new_cond1_final_cc.columns:
        cond1_final_cc_column = new_cond1_final_cc[column].dropna()
        cond2_final_cc_column = cond2_final_cc[column].dropna()
        t_stat, t_p_value = stats.ttest_ind(cond1_final_cc_column, cond2_final_cc_column, equal_var=False, nan_policy='omit')
        new_value = {f"{column} of {condition_1} and {condition_2}": t_p_value}
        new_value_df = pd.DataFrame.from_dict(new_value, orient='index', columns=['p-value'])
        test_values = pd.concat([test_values, new_value_df], axis=0)
    test_values.index.name = 'Significance comparison'
    
    # compute for Independet t-test between columns of condition_1 data
    '''The independent t-test uses statistical test of unequal variance'''
    test_own_values_cond1 = pd.DataFrame()
    for i, col1 in enumerate(cond1_final_cc):
        for col2 in cond1_final_cc[i+1:]:
            cond1_final_cc_col1 = cond1_final_cc[col1].dropna()
            cond1_final_cc_col2 = cond1_final_cc[col2].dropna()
            t_stat_1, t_p_value_1 = stats.ttest_ind(cond1_final_cc_col1, cond1_final_cc_col2, equal_var=False, nan_policy='omit')
            own_new_values = {f"{col1} and {col2} of {condition_1}": t_p_value_1}
            own_new_values_df = pd.DataFrame.from_dict(own_new_values, orient='index', columns=['p-value'])
            test_own_values_cond1 = pd.concat([test_own_values_cond1, own_new_values_df], axis=0)
    test_own_values_cond1.index.name = 'Significance comparison'    
    
    # compute for Independet t-test between columns of condition_2 data
    '''The independent t-test uses statistical test of unequal variance'''
    test_own_values_cond2 = pd.DataFrame()
    for i, col1 in enumerate(cond2_final_cc):
        for col2 in cond2_final_cc[i+1:]:
            cond2_final_cc_col1 = cond2_final_cc[col1].dropna()
            cond2_final_cc_col2 = cond2_final_cc[col2].dropna()
            t_stat_2, t_p_value_2 = stats.ttest_ind(cond2_final_cc_col1, cond2_final_cc_col2, equal_var=False, nan_policy='omit')
            own_new_values_cond2 = {f"{col1} and {col2} of {condition_2}": t_p_value_2}
            own_new_values_cond2_df = pd.DataFrame.from_dict(own_new_values_cond2, orient='index', columns=['p-value'])
            test_own_values_cond2 = pd.concat([test_own_values_cond2, own_new_values_cond2_df], axis=0)
    test_own_values_cond2.index.name = 'Significance comparison'
    
    # save dataframes to csv
    test_values.to_csv(os.path.join(folder_output, (f"t-test between {condition_1} and {condition_2} (C_C).csv")))
    test_own_values_cond1.to_csv(os.path.join(folder_output, (f"t-test between embryonic stage of {condition_1} (C_C).csv")))
    test_own_values_cond2.to_csv(os.path.join(folder_output, (f"t-test between embryonic stage of {condition_2} (C_C).csv")))
    
    return None

#######################################################################################################################################################

def ttestPoletoChromosome(folder_input, filedata_input, folder_output, name_cond1_fileinput, name_cond2_fileinput, condition_1, condition_2):
    # list for files with different conditions
    cond1_final_pc = []
    cond2_final_pc = []
    
    for file in filedata_input:
        filepath = os.path.join(folder_input, file)
        filename = os.path.basename(filepath) # return the name of each file in the folder
        
        # read out the file of condition 1
        if filename == name_cond1_fileinput:
            
            # check for the appropriate format of the csv files
            try:
                df = pd.read_csv(filepath, index_col=None, encoding='utf-8')
            except:
                df = pd.read_csv(filepath, index_col=None, encoding='latin-1')
            cond1_final_pc.append(df)
            
        # read out the file of condition 2
        elif filename == name_cond2_fileinput:
            
            # check for the appropriate format of the csv files 
            try:
                df1 = pd.read_csv(filepath, index_col=None, encoding='utf-8')
            except:
                df1 = pd.read_csv(filepath, index_col=None, encoding='latin-1')
            cond2_final_pc.append(df1)

    # create new dataframes
    cond1_final_pc = pd.concat(cond1_final_pc, ignore_index=True)
    cond2_final_pc = pd.concat(cond2_final_pc, ignore_index=True)
    new_cond1_final_pc = cond1_final_pc.iloc[:,1:] # omit the first column '1-cell'
    
    # compute for Independent t-test between the 2 given conditions
    '''The independent t-test uses statistical test of unequal variance'''
    test_values = pd.DataFrame()
    for column in new_cond1_final_pc.columns:
        cond1_final_pc_column = new_cond1_final_pc[column].dropna()
        cond2_final_pc_column = cond2_final_pc[column].dropna()
        t_stat, t_p_value = stats.ttest_ind(cond1_final_pc_column, cond2_final_pc_column, equal_var=False, nan_policy='omit')
        new_value = {f"{column} of {condition_1} and {condition_2}": t_p_value}
        new_value_df = pd.DataFrame.from_dict(new_value, orient='index', columns=['p-value'])
        test_values = pd.concat([test_values, new_value_df], axis=0)
    test_values.index.name = 'Significance comparison'
    
    # compute for Independet t-test between columns of condition_1 data
    '''The independent t-test uses statistical test of unequal variance'''
    test_own_values_cond1 = pd.DataFrame()
    for i, col1 in enumerate(cond1_final_pc):
        for col2 in cond1_final_pc[i+1:]:
            cond1_final_pc_col1 = cond1_final_pc[col1].dropna()
            cond1_final_pc_col2 = cond1_final_pc[col2].dropna()
            t_stat_1, t_p_value_1 = stats.ttest_ind(cond1_final_pc_col1, cond1_final_pc_col2, equal_var=False, nan_policy='omit')
            own_new_values = {f"{col1} and {col2} of {condition_1}": t_p_value_1}
            own_new_values_df = pd.DataFrame.from_dict(own_new_values, orient='index', columns=['p-value'])
            test_own_values_cond1 = pd.concat([test_own_values_cond1, own_new_values_df], axis=0)
    test_own_values_cond1.index.name = 'Significance comparison'    
    
    # compute for Independet t-test between columns of condition_2 data
    '''The independent t-test uses statistical test of unequal variance'''
    test_own_values_cond2 = pd.DataFrame()
    for i, col1 in enumerate(cond2_final_pc):
        for col2 in cond2_final_pc[i+1:]:
            cond2_final_pc_col1 = cond2_final_pc[col1].dropna()
            cond2_final_pc_col2 = cond2_final_pc[col2].dropna()
            t_stat_2, t_p_value_2 = stats.ttest_ind(cond2_final_pc_col1, cond2_final_pc_col2, equal_var=False, nan_policy='omit')
            own_new_values_cond2 = {f"{col1} and {col2} of {condition_2}": t_p_value_2}
            own_new_values_cond2_df = pd.DataFrame.from_dict(own_new_values_cond2, orient='index', columns=['p-value'])
            test_own_values_cond2 = pd.concat([test_own_values_cond2, own_new_values_cond2_df], axis=0)
    test_own_values_cond2.index.name = 'Significance comparison'
    
    # save dataframes to csv
    test_values.to_csv(os.path.join(folder_output, (f"t-test between {condition_1} and {condition_2} (P_C).csv")))
    test_own_values_cond1.to_csv(os.path.join(folder_output, (f"t-test between embryonic stage of {condition_1} (P_C).csv")))
    test_own_values_cond2.to_csv(os.path.join(folder_output, (f"t-test between embryonic stage of {condition_2} (P_C).csv")))
    
    return None
    