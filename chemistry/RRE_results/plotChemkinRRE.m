function [t,temperature]=plotChemkinRRE(folder)

% This MATLAB function is designed to plot the output files of CHEMKIN II
% when it has calculated the evolution of a reactive fluid particle

% Move to selected folder
prevFolder = cd(folder);
% Extract values from output files of CHEMKIN II
[~, M1]=hfile('CV1.dat');
[~, M2]=hfile('Exp.dat');
[~, M3]=hfile('CV2.dat');

% Construction of time array
timeZ1 = M1(:,1);
timeZ2 = M2(:,1)+timeZ1(end);
timeZ3 = M3(:,1)+timeZ2(end);
t = cat(2, cat(2, timeZ1', timeZ2'), timeZ3');

% Construction of temperature array
tempZ1 = M1(:,2);
tempZ2 = M2(:,2);
tempZ3 = M3(:,2);
temperature = cat(2, cat(2, tempZ1', tempZ2'), tempZ3');

cd(prevFolder);
end