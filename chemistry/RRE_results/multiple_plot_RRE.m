% This MATLAB script is designed to plot several particles evolutions or
% several mach numbers evolution for a RRE structure

run H2_O2_conds.m
omega_deg = 14.5;
temp_I = 600; % in K
temp_II = 1138; % in K
press = 101325; % in Pa
cst = [gamma_I, gamma_II, R_I, R_II];

mach_numbers = [2 2.2 2.4 2.6 2.8 3];
%mach_numbers = 2;
mach_flow = mach_numbers./sin(omega_deg*pi/180);
X = [5 10 15 20 25];
Y = round(X*tan(omega_deg*pi/180),2); % Particles are aligned with the incident shock

X = (X+5)*1e-2;
Y = Y*1e-2;

nM = length(mach_numbers);
nP = length(X);

mode = 2;

colors = [0.6 0 0.8; 0.4 0 0.6; 0.5 0 0.4; 0.2 0.2 0.5; 0 0.4 0.8; 0 0.6 0.5; 0 0 1];

%% Multiple particles
if mode == 1 
    mach = 1;
    
    for part = 1:1:nP
        [sol, expan] = RREevol(cst, X(part), Y(part), mach_numbers(mach), ...
            omega_deg, press, temp_I, temp_II, 1, colors(part+2,:));
    end

%% Multiple mach numbers
else
    part = 1;
    
    for mach = 1:1:nM
        [sol, expan] = RREevol(cst, X(part), Y(part), mach_numbers(mach), ...
            omega_deg, press, temp_I, temp_II, 1, colors(mach+1,:));
    end
end
