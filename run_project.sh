#!/bin/bash

# Compile the C code
gcc -O2 -o tumor_simulation C_code.c -lm

# Run the C simulation
./tumor_simulation

# Run the MATLAB script (without opening GUI)
matlab -nodisplay -nosplash -r "run('Matlab_code.m'); exit;"
