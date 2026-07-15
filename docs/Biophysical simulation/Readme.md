# Biophysical simulation for spindle positioning and elongation

This repository contains a biophysical simulation script used to model and analyze spindle elongation and centrosome positioning dynamics, specifically corresponding to **Figure 6, S7, and S8** of the associated publication.

---

## Description

This script implements a 3D biophysical simulation of mitotic spindle positioning and elongation (typically studied in early embryos like *C. elegans*). It models the cell boundary as a superellipsoid surface reconstructed from a triangulated mesh (`unit_sphere.mat`). 

### Core Biophysical Mechanics
1. **Cell Surface Representation**:
   - The cell boundary is loaded from a triangulated surface (`unit_sphere.mat`), scaled by specified cell radii, and projected onto a superellipsoid geometry.
   - For each surface element (triangle), the incenter, area, and surface normal vector are calculated.

2. **Microtubule Impingement**:
   - Microtubules nucleate from two centrosomes (spindle poles) and grow toward the cell cortex.
   - The rate of microtubule impingement on each cortical surface element is computed based on distance, microtubule growth velocity (`mt_vg`), average length (`mt_ave_length`), nucleation rate (`mt_nuc_rate`), and boundary orientation relative to the centrosome.

3. **Cortical Force-Generators (CFGs)**:
   - Pulling force-generators are distributed across the cell cortex.
   - Microtubules dynamic binding/attachment probability to CFGs is tracked for each surface element. The attachment rate depends on microtubule impingement, and detachment is governed by a constant rate (`cfg_k`).

4. **Force and Motion Dynamics**:
   - Cortical force-generators exert pulling forces on the centrosomes along the direction of the interacting microtubules.
   - The spindle elongates according to an exponential model representing internal driving forces $`l_{dot} = chrom_a (\frac{exp^{(-t / chrom_T)}}{chrom_T})`$.
   - The motion of the two centrosomes is determined by balancing the net pulling forces against viscous drag and coupling interactions (governed by translation viscosity `eta` and elongation resistance `nu`).

5. **Parameter Optimization Support**:
   - The script provides helper functions (`simulate_simulation` and `simulation_error`) designed to interface with optimization algorithms (e.g., `scipy.optimize.differential_evolution`) to fit the biophysical parameters to experimental spindle length time-series data.

---

## Key Task

When executed, the script performs the following sequential tasks:
1. **Initialize Parameters**: Sets up physical constants and biological parameters (microtubule dynamics, cortical force-generator properties, embryo dimensions, viscosity coefficients, and simulation time steps).
2. **Load and Process Geometry**: Loads the raw triangulated mesh from `unit_sphere.mat`, scales it to the target cell size, computes geometric centers, normals, and areas for all triangles, and projects them to form a superellipsoid.
3. **Equilibrate Binding Probabilities**: Runs a pre-equilibration loop (10,000 steps) to determine steady-state microtubule-cortex attachment probabilities prior to starting the dynamic movement simulation.
4. **Run Simulation Loop**: Simulates the system over time (from \(t = 0\) to \(t = 206\) seconds with \($\Delta$ t = 0.1\) seconds), updating at each step:
   - Spindle elongation rate.
   - Geometric vectors and distances between centrosomes and surface elements.
   - Pulling forces acting on both centrosomes.
   - Centrosome velocities and positions (using Euler integration).
   - Microtubule impingement and cortical binding probabilities.
5. **Compare with Experimental Data**: Loads experimental measurements of spindle length over time from `2cell.txt`.
6. **Generate Visualizations and Save Files**: Plots the simulated trajectory against the experimental data and exports data and 3D surface files.

---

## Output

Running the script produces the following outputs in the working directory:

| Output File | Format | Description |
| :--- | :--- | :--- |
| `simulation_results_2cell.png` | Image (PNG) | A 1x2 panel plot containing:<br>1. **Centrosome position *vs.* Time** (independent trajectories of both spindle poles).<br>2. **Spindle Length *vs.* Time** (comparison of simulated spindle length vs. experimental data points from `ncell.txt`). |
| `spindle_length_2cell.csv` | CSV Data | A time-series spreadsheet containing columns: `Time [s]` and `Spindle Length [micron]`. |
| `binding_2cell_t200.ply` | 3D Mesh (PLY) | A polygon file format representing the cell surface mesh. The vertices are color-coded based on the average microtubule binding probability at \(t = 200\text{ s}\) (within a normalized range of 0.3 to 0.7), allowing visualization of cortical force-generator recruitment patterns in 3D rendering software (e.g., MeshLab or PyMOL). |
