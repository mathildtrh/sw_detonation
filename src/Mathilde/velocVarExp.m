function dun = velocVarExp (un, M2, dM2, gamma)
% Calculates the local variation of normal velocity of flow inside the expansion
% knowing that the evolution is isentropic and knowing heat ratio gamma, 
% current velocity u, current Mach number M and Mach number variation dM

dun = un/2 / (1 + (gamma-1)/2 * M2) * dM2;
end