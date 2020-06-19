function pressJump = machtoPressJump(M, gamma, phi)
% Calculates the pressure jump through an oblique jump knowing mach M, heat
% ratio gamma, angle of incidence phi
normalMach = M*sin(phi);
pressJump = (1-gamma+2*gamma*normalMach^2)/(gamma+1);
end