#!/bin/bash

# Supprimer d'anciens fichiers pour repartir propre
rm -f healthy_results.csv smoker_results.csv tumor_evolution.gif
rm -rf frames_lungs

# Compile the C code
gcc -O2 -o tumor_simulation C_code.c -lm

# Run the C simulation
./tumor_simulation

# Run the MATLAB script (without opening GUI)
matlab -nodisplay -nosplash -r "run('Matlab_code.m'); exit;"
