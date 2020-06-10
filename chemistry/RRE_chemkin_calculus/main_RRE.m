% This programm executes series of instructions in order to compute
% chemistry calculus for different lagragian particles 
% and different mach numbers

% It is inspired from Yann's programm for lagrangian particle evolution
% (corrected for specific volume calculation)
% and from Remy's routine to make chemistry calculus with SENKIN from
% CHEMKIN II package

%% Initialisation
% WARNING : here is the absolute path of the directory that contains bash
% files and matlab files : should be changed is machine is changed

cd '/home/mathilde/Documents/PRe_Detonation/sw_detonation/chemistry/RRE_chemkin_calculus';

%% Definition of Lagragian particles and mach numbers

mode = 1; %H2_O2 - He interface
XH2=0.6667;
XO2=0.3333;

mach_numbers = 3;

% X = [-0.9, 0.5, 1, 1.5, -0.5];
% Y = [3, 3, 7, 11, 12];
% X = [0.5, 1, 1.5];
% Y = [3, 5, 8];

X = 0.5;
Y = 3;

nM = length(mach_numbers);
nP = length(X);

%% Directories' names

for m = 1:nM
    for p = 1:nP
        name = strcat("mach", num2str(mach_numbers(m)),...
            "_x", num2str(X(p)), "_y", num2str(Y(p)));
        dirnames(nP*(m-1)+p) = strrep(name,'.',[]); %rises a warning : doesn't matter
    end
end

X=X*1e-2;
Y=Y*1e-2;

%% Main loop for mach numbers
for mach=1:nM
%% Second loop for Lag particles
    for part = 1:nP
        dir = dirnames(nP*(mach-1)+part);
        % change directory to the right one
        mkdir(dir);
        prevDir = cd(dir);
        newfolder = strcat('/home/mathilde/Documents/PRe_Detonation/sw_detonation/chemistry/RRE_chemkin_calculus','/', dir);
        addpath(newfolder);
        % computation of useful data
            % tfinal1 = time at which expansion fan reaches particle
            % tfinal2 = time at which particle quits expansion fan
            % tfinal3 = end of computation
            % P = pressure after incident shock
            % T = temperature after incident shock
        [t, V, tfinal1, tfinal2, tfinal3, P, T] = partThroughRRE(mode, X(part), Y(part),mach_numbers(mach));
        % write VTIM.dat file in the right dir
        file = fopen('VTIM.dat','w');
        fprintf(file, '%1.4E   %1.5f\n', cat(1,t,V));
        fclose(file);

        % Make input and calculation for first constant volume calculations
        
        run makesenkCV1.m
        command = strcat('../bash_CV1.sh', " ", newfolder);
        system(command);

        % Make input and calculation for Expansion calculations
        run makesenkExp.m
        command = strcat('../bash_Exp.sh', " ", newfolder);
        system(command);

        % Make input and calculation for second constant volume calculations
        run makesenkCV2.m
        command = strcat('../bash_CV2.sh'," ", newfolder);
        system(command);
        
        % go back to the main directory
        cd(prevDir)
    end
end