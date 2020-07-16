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
mode = 1; %H2_O2 - He interface

if mode == 1
    run H2_O2_conds.m
    %% Definition of Lagragian particles and mach numbers
    omega_deg = 21.5;
    temp_I = 600; % in K
    temp_II = 1138; % in K
    press = 101325; % in Pa
    cst = [gamma_I, gamma_II, R_I, R_II]; 
    
    % Choose range of pressure jumps to explore
    press_jumps = linspace(1.02, 1.66, 33);
    mach_numbers = sqrt(xiToSqMach(press_jumps, gamma_I, pi/2));
    mach_flow = sqrt(xiToSqMach(press_jumps, gamma_I, omega_deg*pi/180));
    
%     % Choose range of Mach numbers of the incident flow
%     mach_numbers = linspace(3.66, 4.02, 37);
%     mach_flow = mach_numbers./sin(omega_deg*pi/180);
    
    X = [1 5 10 15 20 25];
    Y = round(X*tan(omega_deg*pi/180),2); % Particles are aligned with the incident shock
    X = X+5;
    
    % Choose particles to study
    nM = length(mach_numbers);
    nP = length(X);
    
    %% Directories' names
    
    for m = 1:nM
        for p = 1:nP
            name = strcat("mach", num2str(mach_flow(m)),...
                "_x", num2str(X(p)));
            dirnames(nP*(m-1)+p) = strrep(name,'.','-'); %rises a warning : doesn't matter
        end
    end
    
    X=X*1e-2;
    Y=Y*1e-2;
    
    %% Section to run only for chemical calculus
    if exist("chemical_run","var") ~= 1
        chemical_run = true;
    end
    
    if chemical_run
        %% Main loop for mach numbers
        for mach=1:nM
            %% Second loop for Lag particles
            for part = 1:nP
                out = ['mach = ', num2str(mach_flow(mach)), '; particle x = ', num2str(X(part)), ', y = ', num2str(Y(part))];
                dir = dirnames(nP*(mach-1)+part);
                % change directory to the right one
                mkdir(dir);
                prevDir = cd(dir);
                newfolder = strcat('/home/mathilde/Documents/PRe_Detonation/sw_detonation/chemistry/RRR_chemkin_calculus','/', dir);
                addpath(newfolder);
                % computation of useful data
                % sol = solution of computation for all zones
                sol=RRRevol(cst, X(part), Y(part), mach_numbers(mach), omega_deg, press, temp_I, temp_II, 0);
                
                if sol == -1
                    disp("Solution doesn't exist. EXIT.")
                    cd(prevDir)
                    disp(out);
                    continue
                end
                % Reinitialisation of molar fractions
                XH2=2/3;
                XO2=1/3;
                
                P2 = sol(2,2)/ 101325;
                T2 = sol(3,2);
                P3 = sol(2,3)/ 101325;
                T3 = sol(3,3);
                t_2 = sol(1,2);
                t_3 = sol(1,3);
                % Make input and calculation for first constant volume
                % calculations (after incident shock, zone 2)
                run makesenkzone2.m
                command = strcat('../bash_zone2.sh', " ", newfolder);
                system(command);
                
                % Make input and calculation for second constant volume
                % calculations (after reflected shock, zone 3)
                run makesenkzone3.m
                command = strcat('../bash_zone3.sh'," ", newfolder);
                system(command);
                
                % go back to the main directory
                cd(prevDir)
                disp(out);
            end
        end
    end
end
