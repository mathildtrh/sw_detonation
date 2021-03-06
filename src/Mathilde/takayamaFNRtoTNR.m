% Calculates the angle where FNR-TNR occurs
% ... according to Zeng and Takayama, 1996
function [omega_crit, Msj,Mst,Msi]= takayamaFNRtoTNR(xi_i,angle_gam,gamma_I,...
    gamma_II,mu_I,mu_II,varargin)

%Inputs : 

%xi_i : incident pressure jump
% angle_gam : angle between LIp and horizontal line
% gammas : ratios of specific heats
% mus : molecular masses
T0 = 298;
temp_ratio = 1;
if nargin > 7
    temp_ratio = varargin(1);
end


R = 8.314; %universal constant for perfect gases
a_I = sqrt(gamma_I*R/mu_I*T0); %speed of sound in phase I
a_II = sqrt(gamma_II*R/mu_II*T0*temp_ratio); %speed of sound in phase II
speed_ratio = a_II/a_I;

M1n=sqrt(xiToSqMach(xi_i,gamma_I,pi/2)); %incident shock Mach
M1=M1n/sin(angle_gam); %incident free stream Mach

Msi = M1/sqrt(1+(tan(pi-angle_gam)));

c=(gamma_II+1)/(gamma_I+1)*sqrt(temp_ratio*(gamma_I*mu_II)...
    /(gamma_II*mu_I))*(Msi^2-1)/Msi;
Mst=1/2*(c+sqrt(c^2+4)); %transmitted shock Mach

Vs = a_I*Msi; %speed of the shock
Vt = Mst*a_II;

% Msj is the solution of equation (21) in the paper
% has to be solved like in isBefore FPRToTNR

syms a
eqn= 2*gamma_II/(gamma_I - 1)*Mst^2 - (gamma_II - 1)/(gamma_II + 1) ...
    == (2*gamma_I*a - (gamma_I - 1))/(gamma_I + 1)...
    *(1 - (gamma_II - 1)/(gamma_I + 1)*1/speed_ratio*(a-1/a))^(-2*gamma_II/(gamma_II-1));
res=[(vpasolve(eqn,a,[1,Inf]))];
Msj = res(1);

bet = asin(Msj/Mst*1/speed_ratio); % angle between j-shock and pre-shock flow


omega_crit = pi/2 - angle_gam + asin(Vt*sin(angle_gam)*sin(bet)/Vs) - bet;
end
