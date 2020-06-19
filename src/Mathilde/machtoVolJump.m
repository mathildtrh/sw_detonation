function volJump = machtoVolJump(M, gamma, phi)
% Calculates the specific volume jump through an oblique jump knowing mach M, heat
% ratio gamma, angle of incidence phi
normalMach = M*sin(phi);
volJump = ((gamma+1)*normalMach^2)/(2+(gamma-1)*normalMach^2);
end