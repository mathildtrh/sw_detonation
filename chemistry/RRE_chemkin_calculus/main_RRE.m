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
% Run initial conditions
if mode == 1
    run H2_O2_conds.m
    
    omega_deg = 18;
    temp_I = 600; % in K
    temp_II = 1138; % in K
    press = 101325; % in Pa
    cst = [gamma_I, gamma_II, R_I, R_II];
    % mach_numbers = linspace(3.3, 6.5,20);
    % X = [2, 5.5, 7.5, 10, 12, 13.5];
    % Y = [0.1, 1.78, 0.2, 3, 0.5, 4.38];
    
    mach_numbers = 8;
    X = 9;
    Y = 2;
    
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
            [sol, expan] = RREevol(cst, X(part), Y(part), mach_numbers(mach), omega_deg, press, temp_I, temp_II,0);
            t = expan(1,:);
            V = expan(4,:);
            % write VTIM.dat file in the right dir
            file = fopen('VTIM.dat','w');
            fprintf(file, '%1.4E   %1.5f\n', cat(1,t,V));
            fclose(file);
            
            P = sol(2,2)/101325;
            T = sol(3,2);
            tfinal1 = sol(1,2);
            tfinal2 = t(end);
            tfinal3 = sol(1,3);
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
            
            out = ['mach = ', num2str(mach_numbers(mach)), '; particle x = ', num2str(X(part)), ', y = ', num2str(Y(part))];
            disp(out);
        end
    end
    
else
    disp ("Unexisting mode. EXIT");
end
