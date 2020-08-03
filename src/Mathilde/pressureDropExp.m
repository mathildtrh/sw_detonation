function press_ratio = pressureDropExp(M1, M2, gamma)
% computes the pressure drop through a prandtl meyer expansion

temp_ratio = (1 + (gamma-1)/2*M1^2) / (1 + (gamma-1)/2*M2^2);
press_ratio = temp_ratio^( gamma/(gamma-1) );
end