function dtheta = deflectVarExp (M2, u, du)
% Calculates the local deflection angle inside the expansion
% knowing that the evolution is isentropic 
% and knowing current Mach number M, current velocity of the flow u and
% local variation of velocity

dtheta = sqrt(M2 - 1)*du/u;
end