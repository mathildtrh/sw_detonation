function M1n = machtoPSNormalMach(M, gamma, phi)
% Calculates the post shock normal mach number through an oblique jump knowing mach M, heat
% ratio gamma, angle of incidence phi

% If normal Mach is already known, set M = normalMach and phi = pi/2
normalMach = M*sin(phi);
M1n = sqrt( (1+(gamma+1)/2*normalMach^2) / (gamma*normalMach^2 -(gamma-1)/2) );
end