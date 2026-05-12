# Statistical Analysis

This folder contains Python scripts and Jupyter notebooks for performing statistical validation of the spindle dynamics data.

## 🛠️ Tools & Scripts

### 1. T-Tests
- **`ttest_statistics.ipynb`**: Interactive notebook for calculating independent t-tests on distance analysis data (Pole-to-Pole, Pole-to-Chromosome, Chromosome-to-Chromosome).
- **`testStatistics.py`**: A comprehensive library of t-test functions designed for populations with unequal variance. It compares data between different embryonic stages and conditions (Wild-type vs. RNAi).
- **`welch_ttest_analysis.py`**: Specifically implements Welch's t-test for comparing spindle lengths across different cell types, providing detailed pairwise comparison results.

### 2. ANOVA & Mixed Linear Models
- **`ANOVA_Mixed_Linear_Model.ipynb`**: Notebook for advanced statistical modeling using Mixed Linear Models to account for variance components across stages and cells.
- **`anova.py`**: Python script using `statsmodels` to perform Mixed Linear Model analysis on spindle length data, extracting variance components for measurement noise and biological variability.

## 📊 Input Data Structure
Most scripts expect CSV files with columns like `Stage`, `Cell`, and the measurement of interest (e.g., `Final pole_pole length (µm)`).
For t-tests, the input should often be structured with stages as columns:
| 1-cell | 2-cell | 4-cell | ... |
| :---: | :---: | :---: | :---: |
| val | val | val | ... |

## 🚀 How to Run
1. Ensure `scipy`, `statsmodels`, `pandas`, and `numpy` are installed.
2. For notebooks:
   ```bash
   jupyter notebook ttest_statistics.ipynb
   ```
3. For scripts:
   ```bash
   python welch_ttest_analysis.py
   ```
   (Note: You may need to edit the script to point to your specific data file).
