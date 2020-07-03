function pressJump = machtoPressJump(M, gamma, phi)
% Calculates the pressure jump through an oblique jump knowing mach M, heat
% ratio gamma, angle of incidence phi

% If normal Mach is already known, set M = normalMach and phi = pi/2
normalMach = M*sin(phi);
pressJump = (1-gamma+2*gamma*normalMach^2)/(gamma+1);
end