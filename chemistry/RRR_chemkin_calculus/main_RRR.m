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

cd '/home/mathilde/Documents/PRe_Detonation/sw_detonation/chemistry/RRR_chemkin_calculus';

%% Definition of Lagragian particles and mach numbers

mode = 1; %H2_O2 - He interface
XH2=0.6667;
XO2=0.3333;

%mach_numbers = linspace(2,6);
mach_numbers = [3 3.2];
X = [0.5, 1.5, 0.5, 3, 0.5, 4];
Y = [2, 5.5, 7.5, 10, 12, 13.5];

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
        newfolder = strcat('/home/mathilde/Documents/PRe_Detonation/sw_detonation/chemistry/RRR_chemkin_calculus','/', dir);
        addpath(newfolder);
        % computation of useful data
            % tfinal1 = time at which expansion fan reaches particle
            % tfinal2 = time at which particle quits expansion fan
            % tfinal3 = end of computation
            % P = pressure after incident shock
            % T = temperature after incident shock
        [P0, T0, t_0, P1, T1, t_1, P2, T2, t_2] = RRRevolution(mode, X(part), Y(part),mach_numbers(mach),false);

        % Make input and calculation for first constant volume calculations
        
        run makesenkzone0.m
        command = strcat('../bash_zone0.sh', " ", newfolder);
        system(command);

        % Make input and calculation for second constant volume calculations
        run makesenkzone1.m
        command = strcat('../bash_zone1.sh', " ", newfolder);
        system(command);

        % Make input and calculation for third constant volume calculations
        run makesenkzone2.m
        command = strcat('../bash_zone2.sh'," ", newfolder);
        system(command);
        
        % go back to the main directory
        cd(prevDir)
        
        out = ['mach = ', num2str(mach_numbers(mach)), '; particle x = ', num2str(X(part)), ', y = ', num2str(Y(part))];
        disp(out);
    end
end
