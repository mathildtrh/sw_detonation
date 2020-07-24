% This programm computes the evolution of pressure, temperature and
% specific volume of a Lagrangian particle when it travels through a RRE
% structure

% Inputs : mode = 1 : H2_O2 // He interface
%          x, y coordinates of the lag. particle at the beginning of the
%          simulation
%          mach : mach number of the incident shock
%          plot : if true, plots the evolution of each variable
% Outputs : shock : true if the reflected wave is a shock (P3>P2)

function [shock]=BPRdisc(cst, mach, omega_deg, press, temp_I, temp_II)

shock = false;

gamma_I = cst(1);
gamma_II = cst(2);
R_I = cst(3);
R_II = cst(4);

omega_i = omega_deg*pi/180; % angle of incidence in rad

%% Zone 1

P1 = press;                                                            % Pressure
T1 = temp_I;                                                           % Temperature
V1 = R_I*T1/P1;                                                        % Spec. volume
a1 = sqrt(gamma_I*R_I*T1);                                             % Speed of sound

%% Change of frame of reference

ui = mach*a1;
M1 = mach/sin(omega_i);                                                % Mach number in zone 1

%% Zone 2 (oblique shock from 1)

xi_i = machtoPressJump(M1, gamma_I, omega_i);                          % Pressure jump
P2 = xi_i*P1;                                                          % Pressure

%% Zone 5

P5 = press;                                                            % Pressure
T5 = temp_II;                                                          % Temperature
a5 = sqrt(gamma_II*R_II*T5);                                           % Speed of sound

%% Zone 4 (oblique shock from 5)
b = (gamma_II+1)/(gamma_I+1)*(ui^2-a1^2)/ui;                           % Temporary value
ut = 0.5*(b+sqrt(b^2+4*a5^2));                                         % Normal speed of transmitted shock
omega_t = asin(sin(omega_i)*ut/ui);                                    % Angle of transmitted shock
M5n = ut/a5;                                                           % Normal Mach number
M5 = M5n/sin(omega_t);                                                 % Mach number

xi_t = machtoPressJump(M5, gamma_II, omega_t);                          % Pressure jump
P4 = xi_t*P5;                                                          % Pressure

%% Zone 3 (membrane equilibrium with 4)

P3 = P4;                                                               % Pressure P3 = P4 = P5*xi_t
xi_r = P3/P2;                                                          % Reflected pressure jump
if xi_r > 1
    shock = true;
end

end
