function dV = volVarExp (V, M2, dM2, gamma)
% Calculates the local variation of spec. volume inside the expansion
% knowing that the evolution is isentropic and knowing heat ratio gamma, 
% current spec. volume V, current Mach number M and Mach number variation dM

dV = V * dM2/2 / (1 + (gamma-1)/2 * M2);

end