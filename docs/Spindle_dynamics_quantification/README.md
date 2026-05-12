# Spindle Dynamics Quantification & Visualization

This folder contains Jupyter notebooks dedicated to the quantification of spindle elongation and chromosome segregation dynamics, as well as the generation of figures for publication.

## 📓 Notebooks

### 1. Spindle Elongation (Pole-to-Pole)
- **`Figure_1f_mean_P_P_curve_fitting_SIGMOID_function.ipynb`**:
    - **Purpose**: Fits a sigmoid function to average pole-to-pole distance data.
    - **Output**: Extracted parameters (Initial/Final length, Metaphase length, Elongation rate) and Fig. 1f plots.

### 2. Chromosome Segregation (Chromosome-to-Chromosome)
- **`Figure_1g_mean_C_C_curve_fitting_EXPONENTIAL_function.ipynb`**:
    - **Purpose**: Fits an exponential function to average chromosome-to-chromosome distance data.
    - **Output**: Extracted parameters (Final segregation length, Segregation speed) and Fig. 1g plots.

### 3. Dynamics & Comparison Plots
- **`Figure_2c_2g_dynamic_group_plots.ipynb`**: Visualizes average spindle elongation and chromosome segregation over time across different developmental stages and conditions (Wild-type vs. RNAi).
- **`Figure_2d_2e_2g_2i_COMBINED_plots_Spindle_parameters.ipynb`**: Plots final spindle length and chromosome segregation distance as a function of cell size, including linear regression analysis.
- **`Figure_2f_2j_dynamics_summary_plots.ipynb`**: Quantifies and summarizes spindle elongation rates and segregation speeds as a function of cell size across different embryonic conditions.

### 4. Supplemental Analysis
- **`Figure_supplement_6_plot_and_statistics_ELONGATION_RATE.ipynb`**: Detailed analysis and statistical plotting of elongation rates, specifically for supplemental figures.

## 🚀 How to Run
1. Ensure your environment has `pandas`, `seaborn`, `matplotlib`, and `scipy`.
2. Open the notebook corresponding to the figure or analysis you wish to reproduce.
3. Ensure the input CSV files (e.g., `pole_pole_distance.csv`, `chromosome_distance.csv`) are available in the expected data directories.
4. Run the cells to generate the plots and extract statistical parameters.
