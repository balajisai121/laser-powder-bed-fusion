# laser-powder-bed-fusion
A MATLAB-based simulation of thermal fields and melt pool dynamics during Laser Powder Bed Fusion (L-PBF) of Inconel 625. Includes transient thermal analysis, Gaussian heat source modeling, and experimental validation.

# Laser Powder Bed Fusion Simulation

This repository contains a MATLAB project focused on simulating the thermal fields and melt pool dynamics during the Laser Powder Bed Fusion (L-PBF) process for Nickel Alloy 625. The project models heat transfer, material properties, and boundary conditions, providing insights into melt pool dimensions and temperature distributions.

## Project Objectives
- Simulate the **transient thermal field** in Laser Powder Bed Fusion.
- Analyze the **melt pool size, shape, and depth** under varying laser power and scanning speeds.
- Validate results through experimental data.

## Key Features
- **3D Thermal Simulation**: Models heat conduction using a transient heat equation.
- **Gaussian Beam Heat Source**: Realistic energy deposition modeled via a moving Gaussian heat source.
- **Material Properties**: Temperature-dependent thermal properties of Inconel 625.
- **Melt Pool Dynamics**: Analysis of melt pool length, width, and depth for different energy densities.

## Methods
- **Heat Transfer Equation**:
  \[
  \rho C_p \frac{\partial T}{\partial t} - \nabla \cdot (k \nabla T) = Q_v
  \]
  where \( T \) is the temperature, \( \rho \) is density, \( C_p \) is specific heat, \( k \) is thermal conductivity, and \( Q_v \) is the heat source.

- **Gaussian Heat Source**:
  \[
  Q(x, y) = \frac{(1 - R)P}{\pi w_0^2} \exp{\left(-\frac{2r^2}{w_0^2}\right)}
  \]
  where \( P \) is laser power, \( R \) is reflectivity, and \( w_0 \) is beam radius.

- **Boundary Conditions**:
  - **Dirichlet**: Fixed temperature at the base plate (353 K).
  - **Neumann**: Ambient heat flux at exposed surfaces (293 K).

## Simulation Results
- **Melt Pool Dimensions**:
  - **Low Energy Density**: Length: 0.299 mm, Width: 0.164 mm, Depth: 0.030 mm
  - **Medium Energy Density**: Length: 0.350 mm, Width: 0.165 mm, Depth: 0.033 mm
  - **High Energy Density**: Length: 0.425 mm, Width: 0.165 mm, Depth: 0.043 mm

- **Temperature Distribution**:
  - Significant thermal gradients promote directional solidification and affect microstructure growth.

## Visualizations
- 3D plots of temperature distributions for varying energy densities.
- Scatter plots showing temperature variation along scanning paths.
- Melt pool visualization under different process conditions.

## Repository Structure
```plaintext
laser-powder-bed-fusion/
├── data/                 # Simulation data and results
│   ├── melt_pool_data.csv
│   ├── temperature_profiles.csv
│   └── experimental_results.csv
├── scripts/              # MATLAB scripts for simulation
│   ├── thermal_model.m
│   ├── heat_source.m
│   └── boundary_conditions.m
├── notebooks/            # Jupyter notebooks for data analysis
│   ├── data_analysis.ipynb
│   ├── melt_pool_analysis.ipynb
│   └── visualization.ipynb
├── visualizations/       # Plots and images
│   ├── 3d_temperature_plot.png
│   ├── melt_pool_plot.png
│   └── validation_results.png
├── docs/                 # Reports and presentations
│   ├── FINAL_REPORT.pdf
│   ├── FINAL_PRESENTATION.pdf
│   └── midterm_exam.pdf
├── README.md             # Project overview
├── LICENSE               # License file
└── requirements.txt      # Python dependencies for Jupyter notebooks
