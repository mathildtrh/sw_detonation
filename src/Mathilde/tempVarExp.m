function dT = tempVarExp (T, M2, dM2, gamma)
% Calculates the local variation of temperature inside the expansion
% knowing that the evolution is isentropic and knowing heat ratio gamma, 
% current temperature T, current Mach number M and Mach number variation dM

dT = -T * (gamma-1)/2*dM2 / (1 + (gamma-1)/2 * M2);

end