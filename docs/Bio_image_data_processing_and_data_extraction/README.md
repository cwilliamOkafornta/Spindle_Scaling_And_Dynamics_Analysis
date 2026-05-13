# Bio-image data processing and data extraction

This folder contains scripts and macros for processing raw microscopy images (typically multi-channel TIFF files) to extract quantitative data about centrosome and chromosome positions over time.

## 🛠️ Tools & Scripts

### 1. ImageJ Macros (`.ijm`)
- **`spot_finder.ijm`**: Used for automated identification and tracking of fluorescent spots (centrosomes) within Fiji/ImageJ.
- **`kymograph script.ijm`**: Generates kymographs from time-lapse movies to visualize the movement of spindle components along a defined path over time.
    - The kymograph script, in addition to the input image data, uses 3D positions of the centrosomes from the **`spot finder.ijm`** to generate kymographs.

### 2. MATLAB Scripts
- **`centrosome_chromosome_intensity.m`**: The core analysis script. It performs Gaussian fitting to precisely locate centrosome positions and utilizes peak detection to identify chromosome positions within the spindle.
    - **Input**: Multi-channel TIFF images.
        - The folder with the Multi-channel TIFF images should be placed in the same folder as the two MATLAB scripts.
    - **Output**: 
        - CSV file containing time-resolved positions (Pole 1, Pole 2, Chromosome 1, Chromosome 2, and number of spots).
        - PNG plot for visual validation of the tracking.
- **`master_file.m`**: A batch processing wrapper for `centrosome_chromosome_intensity.m`. It iterates through all `.tif` files in a specified directory using parallel processing (`parfor`) to accelerate data extraction.


## 🚀 How to Run

### ImageJ Macros
1. Open Fiji/ImageJ.
2. Drag and drop the `.ijm` script into Fiji.
3. Select your input files or directories and click **Run**.

### MATLAB Scripts
1. Open MATLAB.
2. Ensure the `tiffreadVolume` function is in your path (or install the appropriate toolbox).
3. Run the batch processor by calling:
   ```matlab
   Master('path/to/your/image/folder');
   ```
4. Results will be saved as CSV and PNG files in the same folder as the input images.
