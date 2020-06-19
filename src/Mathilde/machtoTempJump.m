function tempJump = machtoTempJump(M, gamma, phi)
% Calculates the temperature jump through an oblique jump knowing mach M, heat
% ratio gamma, angle of incidence phi
normalMach = M*sin(phi);
tempJump = (1-gamma+2*gamma*normalMach^2)/(gamma+1) ...
    * ((gamma+1)*normalMach^2)/(2+(gamma-1)*normalMach^2);
end