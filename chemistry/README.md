# Chemistry folder

This folder is aimed at putting together the MATLAB scripts and the CHEMKIN II files needed to compute the evolution of a Lagragian particle through a given refraction pattern.

## "Chemkin_calculus" folders

These folders contain all the input files for CHEMKIN to run properly + the MATLAB & bash scripts needed to have a full computation.
Open the `mainRRE.m` or `mainRRR.m` file, choose the range of particles desired and the range of Mach numbers to study and run the script. It will automatically calculate the history of pressure, temperature and specific volume of each particle for each Mach number and run CHEMKIN II to simulate chemical kinetics.
The result is organised by subfolders, one for each particle at each Mach number.

## "Results" folders

These folders contain scripts written to help graphical representation of the results.
Scripts can be run before or after `mainRRE.m` or `mainRRR.m` file, just select the right option.

## The special case of BPR

It was impossible to have relevant results with this refraction structure as it is irregular. The present MATLAB scripts are only here to help understand what is going on in this case.

