#!/bin/bash

# Compile the C program located in src/
gcc -O2 -o src/tumor_simulation src/main.c -lm

# Run the C simulation (output will appear in the project folder)
./src/tumor_simulation

# Run the MATLAB script located in src/
/usr/local/bin/matlab-2021b -batch "run('src/make_gif.m')"
