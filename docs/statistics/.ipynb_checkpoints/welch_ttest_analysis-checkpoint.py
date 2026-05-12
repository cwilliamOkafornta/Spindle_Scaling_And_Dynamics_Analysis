
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from itertools import combinations

def welch_ttest_analysis(data_file_path):
    """
    Perform Welch t-tests on spindle length data

    Parameters:
    data_file_path: str, path to CSV file with columns 'Stage', 'Cell', 'Final pole_pole length (µm)'

    Returns:
    DataFrame with all pairwise comparison results
    """

    # Load data
    df = pd.read_csv(data_file_path)

    # Get unique cell types
    cell_types = df['Cell'].unique()

    # Create results list
    results = []

    # Generate all pairwise combinations
    for cell1, cell2 in combinations(cell_types, 2):
        # Extract data for each cell type
        data1 = df[df['Cell'] == cell1]['Final pole_pole length (µm)'].values
        data2 = df[df['Cell'] == cell2]['Final pole_pole length (µm)'].values

        # Perform Welch t-test (equal_var=False)
        t_stat, p_value = ttest_ind(data1, data2, equal_var=False)

        # Calculate descriptive statistics
        mean1 = np.mean(data1)
        mean2 = np.mean(data2)
        std1 = np.std(data1, ddof=1)
        std2 = np.std(data2, ddof=1)
        n1 = len(data1)
        n2 = len(data2)

        results.append({
            'Cell_Type_1': cell1,
            'Cell_Type_2': cell2,
            'Mean_1': mean1,
            'Std_1': std1,
            'n_1': n1,
            'Mean_2': mean2,
            'Std_2': std2,
            'n_2': n2,
            'T_Statistic': t_stat,
            'P_Value': p_value,
            'Significant_0.05': p_value < 0.05,
            'Significant_0.01': p_value < 0.01
        })

    return pd.DataFrame(results)

# Example usage:
# results = welch_ttest_analysis('your_data.csv')
# significant_results = results[results['P_Value'] < 0.05].sort_values('P_Value')
# print(significant_results)
