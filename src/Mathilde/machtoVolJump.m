function volJump = machtoVolJump(M, gamma, phi)
% Calculates the specific volume jump through an oblique jump knowing mach M, heat
% ratio gamma, angle of incidence phi

% If normal Mach is already known, set M = normalMach and phi = pi/2
normalMach = M*sin(phi);
volJump = (2+(gamma-1)*normalMach^2)/((gamma+1)*normalMach^2);
end