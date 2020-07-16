function [t,temperature]=plotChemkinRRR(folder)

% This MATLAB function is designed to plot the output files of CHEMKIN II
% when it has calculated the evolution of a reactive fluid particle

% Move to selected folder
prevFolder = cd(folder);
% Extract values from output files of CHEMKIN II
[~, M2]=hfile('ZN2.dat');
[~, M3]=hfile('ZN3.dat');

% Construction of time array
timeZ2 = M2(:,1);
timeZ3 = M3(:,1)+timeZ2(end);

t = cat(2, timeZ2', timeZ3');

% Construction of temperature array
tempZ2 = M2(:,2);
tempZ3 = M3(:,2);

temperature = cat(2, tempZ2', tempZ3');

cd(prevFolder);
end