# Fit Function Analysis

This folder contains Jupyter notebooks designed to fit mathematical models to experimental distance data extracted from *C. elegans* embryos. The goal is to determine key kinetic parameters such as initial/final lengths and rates of movement.

## 📓 Notebooks

### 1. Pole-to-Pole Distance Fitting
- **`curve_fit_P_P_dropna_Batch_process(Individual_fit).ipynb`**:
    - **Model**: Sigmoid function: $L = a + \frac{b}{1 + e^{-(t-t0)/\tau}}$
    - **Extracted Parameters**: Initial pole-to-pole length, Final pole-to-pole length, Elongation rate, and Metaphase length.
    - **Usage**: Fits individual experimental datasets (n-values) in batch.

### 2. Chromosome-to-Chromosome Distance Fitting
- **`curve_fit_C_C_dropna_Batch_process(Individual_fit).ipynb`**:
    - **Model**: Exponential function: $y = a(1 - e^{-x/b})$
    - **Extracted Parameters**: Final chromosome-to-chromosome length and Segregation speed.
    - **Usage**: Batch processes individual experiment data for different cell types.

### 3. Pole-to-Chromosome Distance Fitting
- **`curve_fit_Avg_P_C_dropna_Batch_process(Individual_n_fit).ipynb`**:
    - **Model**: Polynomial function: $y = p_n x^n + ... + p_0$
    - **Extracted Parameters**: Initial length, Final length, Rate of distance reduction, and Normalized ratio.
    - **Usage**: Determines how the distance between the pole and chromosomes changes during segregation.

## 📊 Input Data Structure
The notebooks expect `.csv` files with the following structure:
| Exp00 | Exp01 | ... | mean | std | n | SE | time |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| num | num | ... | num | num | num | num | num |

## 🚀 How to Run
1. Ensure you have a Python environment with `numpy`, `pandas`, `matplotlib`, and `scipy` installed.
2. Open the desired notebook in Jupyter Lab or Jupyter Notebook.
3. Update the `folder_path` or `file_path` variables in the first few cells to point to your CSV data.
4. Run all cells to perform the batch fitting and visualize the results.
