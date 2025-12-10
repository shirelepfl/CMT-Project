@echo off
REM Supprimer anciens fichiers pour repartir propre
del /Q healthy_results.csv
del /Q smoker_results.csv
del /Q tumor_evolution.gif
rmdir /S /Q frames_lungs

REM Compiler le code C
gcc -O2 -o tumor_simulation C_code.c -lm

REM Lancer la simulation
tumor_simulation.exe

REM Lancer MATLAB
matlab -nodisplay -nosplash -r "run('Matlab_code.m'); exit;"
