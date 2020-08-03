function temp_ratio = temperatureDropExp(M1, M2, gamma)
% computes the temperature drop through a prandtl meyer expansion
temp_ratio = (1 + (gamma-1)/2*M1^2) / (1 + (gamma-1)/2*M2^2);

end