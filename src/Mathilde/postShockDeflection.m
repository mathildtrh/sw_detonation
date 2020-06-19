function delta = postShockDeflection(M, gamma, phi)
% Calculates the temperature jump through an oblique jump knowing mach M, heat
% ratio gamma, angle of incidence phi
normalMach = M*sin(phi);
tgt = 2*cot(phi)*(normalMach^2 - 1)/(M^2 * (gamma + cos(2*phi)) +2) ;
delta = atan(tgt);
end