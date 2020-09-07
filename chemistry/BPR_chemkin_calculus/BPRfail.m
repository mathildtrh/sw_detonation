%Computation of the evolution of Lagragian particle of coordinates x and y
%in the coordinate system set in theory (see report)
%for a BPR refraction structure

%Stany Gallier's simulation are used to compare results:
%it is found that the incident shock has the expected strength
%but our computation of the reflected shock is wrong if the angle of the
%reflected shock the one measured on Gallier's simulations


function [rho1, rho2a, rho2b] = BPRfail (press)

gamma_I=1.6696;
gamma_II=1.6696;
mu_I=40;
mu_II=4;
Ru = 8.314;
R_I = Ru/(mu_I*1e-3);
R_II = Ru/(mu_II*1e-3);

omega_deg = 27;
omega_i = omega_deg*pi/180; % angle of incidence in rad
chi = 0.131;
xi_i = 1/chi;

%% Zone 1

P1 = press;                                                            % Pressure
rho1 = 1.6145;                                                           % Density
V1 = 1/rho1;                                                           % Spec. volume
M1 = sqrt(xiToSqMach(xi_i, gamma_I, omega_i));                         % Mach number in zone 1

%% Zone 2a (oblique shock from 1)

M2an = machtoPSNormalMach(M1, gamma_I, omega_i);                       % Mach number (2) normal to incident shock
V2a = V1*machtoVolJump(M1, gamma_I, omega_i);                          % Spec. volume
rho2a = 1/V2a;                                                         % Density
delta_1 = postShockDeflection(M1, gamma_I, omega_i);                   % Deflection angle behind the shock
M2a = M2an / sin(omega_i - delta_1);                                   % Mach number

%% Zone 2b is determined by measuring the angle alpha_r between the initial
% interface and the reflected shock on Gallier's simulations

alpha_r = 177*pi/180; 
omega_r = alpha_r - delta_1;

V2b = V2a*machtoVolJump(M2a, gamma_I, pi-omega_r);                     % Spec Volume
rho2b = 1/V2b;                                                         % Density

end
