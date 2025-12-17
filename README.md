# Lung Cancer Modelisation

## Project Description

This project aims to create a visual awareness campaign on smoking prevention for a young audience. We developed a computational model comparing two scenarios: two young individuals of similar age, one athletic and healthy, and the other a smoker with unhealthy habits.

A C program computes cancer cell progression using the KPP (Kolmogorov–Petrovsky–Piskunov) equation, also known as the Fisher equation. A MATLAB script then visualizes the results by modelling the spatial spread of cancer cells in lung tissue.

### Input files

No external datasets were required. All initial conditions and parameter values were based on typical values found in scientific literature and adapted to fit the two lifestyle scenarios.  

### Output files

The C program generates two CSV files, each containing the number of cancer cells at successive time steps (every 5 days over a 120-day period).

The MATLAB script produces a GIF animation showing lung tissue where cells progressively turn red when becoming cancerous. The animations are generated sequentially for both the healthy and smoker scenarios and can be found in MATLAB’s output files (preview mode).

### Report

The report can be found in the docs folder, under the name 'Report.md'. 

## Running the program

### Dependencies

On the C code, the libraries that need to be included are: stdio.h, stdlib.h, math.h and string.h. The version of MatLab we used is the one from 2025. 

### Build
 
The C code is compiled with gcc on Visual Studio Code. To run the C code and the Matlab code, open the src folder from the CMT-Project. 

```
gcc -Wall C_code.c -o project -lm && ./project
```

### Execute

This project consists of a C simulation (C_code.c) that generates CSV data and a MATLAB script (Matlab_code.m) that creates a GIF from the simulation results. The two codes, thus, need to be kept in the same folder. 


1. Open a terminal and clone the repository:


```bash

git clone https://github.com/shirelepfl/CMT-Project.git

cd CMT-Project/src

```

2. Execute the project:

   ```
   make
   ```
The output files can be found in the src folder. 

## Contributors

Shirel Gelfged and Rebecca Levi. 

## Acknowledgments

### Data sources

Initial parameter ranges (net proliferation rate, effective diffusion, carrying capacity, and immune-clearance term) are based on: https://www.researchgate.net/publication/7954500_A_hybrid_mathematical_model_of_solid_tumour_invasion_The_importance_of_cell_adhesion.

Parameter coefficients reflecting lifestyle differences were obtained from:

Smoker lungs:
- the net proliferation rate coefficient value : https://acsjournals.onlinelibrary.wiley.com/doi/10.1002/1097-0142%2820001001%2989%3A7%3C1457%3A%3AAID-CNCR7%3E3.0.CO%3B2-L
- the effective diffusion coefficient value : https://onlinelibrary.wiley.com/doi/10.1155/2019/2025636
- the carrying capacity coefficient value : https://www.sciencedirect.com/science/article/pii/S0024320512003402?via%3Dihub
- the net immune-clearance term coefficient value : https://www.nature.com/articles/s41598-020-76556-7

Healthy lungs:
- the net proliferation rate coefficient value, (same source as smoker lungs) : https://acsjournals.onlinelibrary.wiley.com/doi/10.1002/1097-0142%2820001001%2989%3A7%3C1457%3A%3AAID-CNCR7%3E3.0.CO%3B2-L
- the effective diffusion coefficient value, (same source as smoker lungs) : https://onlinelibrary.wiley.com/doi/10.1155/2019/2025636
- the carrying capacity coefficient value : https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1474770/full
- the net immune-clearance term coefficient value : https://sportsmedicine-open.springeropen.com/articles/10.1186/s40798-022-00419-w


### Code

The code was written with help from the assistants and ChatGPT version 5.1. No pre-existing online codebases were used.
