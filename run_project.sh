#!/bin/bash

# Install essential development libraries (includes stdio, stdlib, math, string)
sudo apt-get update
sudo apt-get install -y build-essential


# Compile the C program located in src/
gcc -O2 -o results/output src/C_code.c -lm

# Run the C simulation (output will appear in the project folder)
./results/output

# Run the MATLAB script located in src/
/usr/local/bin/matlab-2021b -batch "run('src/Matlab_code.m')"
