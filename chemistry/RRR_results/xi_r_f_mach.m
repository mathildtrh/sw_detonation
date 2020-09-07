run H2_O2_conds.m
%% Definition of Lagragian particles and mach numbers
omega_deg = 21.5;
temp_I = 600; % in K
temp_II = 1138; % in K
press = 101325; % in Pa
cst = [gamma_I, gamma_II, R_I, R_II];

% Choose range of pressure jumps to explore
press_jumps = linspace(1.02, 1.66, 23);
mach_numbers = sqrt(xiToSqMach(press_jumps, gamma_I, pi/2));
mach_flow = sqrt(xiToSqMach(press_jumps, gamma_I, omega_deg*pi/180));

nM = length(mach_numbers);

xi_r = zeros(1, nM);
X = 5;
Y = round(X*tan(omega_deg*pi/180),2); % Particles are aligned with the incident shock
X = (X+5)/100;
Y = Y/100;
    
for mach=1:nM
    sol=RRRevol(cst, X, Y, mach_numbers(mach), omega_deg, press, temp_I, temp_II, 0);
    xi_r(mach) = sol(2,3)/sol(2,2);
end

plot(mach_flow,xi_r, 'LineWidth', 2);