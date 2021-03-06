function vol_jump = xiToVolJump(xi, gamma)
% calculates the specific volume jump through the shock ...
% ... knowing pressure jump and specific heat ratio ...
% ... of the phase
    vol_jump = ( (gamma-1)^2 + xi*(gamma^2 -1) + 4*gamma ) / ( gamma^2 - 1 + xi*(gamma+1)^2 );
end