function temp_jump=xiToTempJ(xi,gamma)

% calculates the temperature jump through the shock ...
% ... knowing pressure jump and specific heat ratio ...
% ... of the phase
    temp_jump=xi*((gamma-1)*xi+(gamma+1))/((gamma+1)*xi+(gamma-1));
end