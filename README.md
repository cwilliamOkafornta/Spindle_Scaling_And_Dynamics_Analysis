# *C. elegans* spindle length scaling and chromosome dynamics analysis

This repository contains a comprehensive suite of tools for processing, analyzing, and visualizing spindle length and chromosome dynamics in *C. elegans* embryos. The workflow covers the entire pipeline from raw image data processing to publication-quality statistical analysis.

## 📁 Repository Structure

The core of the analysis pipeline is located in the `docs/` directory, organized into specialized sub-folders:

### 1. [Bio-image Data Processing and Data Extraction](./docs/Bio_image_data_processing_and_data_extraction/)
- **Description**: MATLAB scripts and ImageJ macros for raw TIFF image processing.
- **Key Task**: Extracts centrosome and chromosome positions using Gaussian fitting and peak detection.
- **Output**: Time-resolved distance CSVs and validation plots.

### 2. [Fit Function Analysis](./docs/Fit_function/)
- **Description**: Jupyter notebooks for mathematical modeling of spindle components.
- **Key Task**: Fits Sigmoid, Exponential, and Polynomial curves to individual experiment datasets.
- **Output**: Kinetic parameters like elongation rates, segregation speeds, and initial/final lengths.

### 3. [Spindle Dynamics Quantification](./docs/Spindle_dynamics_quantification/)
- **Description**: Notebooks dedicated to generating publication figures (Figures 1, 2, and supplemental).
- **Key Task**: Quantifies dynamics as a function of cell size and developmental stage.
- **Output**: Publication-quality plots and summary statistics.

### 4. [Microtubule Interactions](./docs/Microtubule_interactions/)
- **Description**: Advanced 3D analysis of microtubule (MT) segments from Amira SpatialGraph data.
- **Key Task**: Calculates interaction proximity, tortuosity, and spatial metrics.
- **Output**: 3D heatmaps, interactive reports, and AmiraMesh exports. Includes a custom **napari plugin**.

### 5. [Statistical Analysis](./docs/statistics/)
- **Description**: Python scripts and notebooks for statistical validation.
- **Key Task**: Performs independent t-tests (Welch's), ANOVA, and Mixed Linear Models.
- **Output**: P-values and variance component analysis for comparing experimental conditions.

## 🚀 Getting Started

### Prerequisites
- **Fiji/ImageJ**: Required for spot-finding and kymograph generation macros.
- **MATLAB**: Required for signal intensity readout from the generated kymographs.
- **Python 3.9+: Recommended to use a virtual environment.
  ```bash
  pip install -r environment.yml  # or use conda
  ```

### General Workflow
1. **Extraction**: Use MATLAB scripts in `Bio_image_data_processing...` to convert `.tif` images to `.csv` data.
2. **Modeling**: Use `Fit_function` notebooks to extract kinetic parameters from the `.csv` files.
3. **Quantification**: Use `Spindle_dynamics_quantification` notebooks to aggregate data and generate figures.
4. **Validation**: Run the scripts in `statistics` to ensure scientific significance of the findings.
5. **Advanced MT Analysis**: For tomography data, use the `Microtubule_interactions` pipeline.

## ⚖️ License
This project is licensed under the BSD-3-Clause License. See the `LICENSE` file for details.
