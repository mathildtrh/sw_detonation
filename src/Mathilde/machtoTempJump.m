function tempJump = machtoTempJump(M, gamma, phi)
% Calculates the temperature jump through an oblique jump knowing mach M, heat
% ratio gamma, angle of incidence phi

% If normal Mach is already known, set M = normalMach and phi = pi/2
normalMach = M*sin(phi);
tempJump = (1-gamma+2*gamma*normalMach^2)/(gamma+1) ...
    * (2+(gamma-1)*normalMach^2)/((gamma+1)*normalMach^2);
end