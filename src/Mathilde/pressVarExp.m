function dP = pressVarExp (P, M2, dM2, gamma)
% Calculates the local variation of pressure inside the expansion
% knowing that the evolution is isentropic and knowing heat ratio gamma, 
% current pressure P, current Mach number M and Mach number variation dM

dP = -P * gamma/2 *dM2 / (1 + (gamma-1)/2 * M2);

end