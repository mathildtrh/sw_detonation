% This program is designed to get the intersection point between two polars
% USEFUL : getPolarPoint.m returns xi=-1 if there is no point coresponding
% to a certain set of coordinates

function [inter_delta, inter_xi] = getCrossPoint(deltas1, xis1, deltas2, xis2)

% while no intersection point is found
% choose a point in the 1st polar (transmitted) = coordinates
% try to find a point in the 2nd polar (reflected) that reaches coordinates
% if there isn't try another point
% if the distance between the meeting point and the coordinate is minimal
% keep this point in memory

len1 = length(deltas1);

n = 1;
coordinate = deltas1(n);
distance = max(xis1); % init of distance between to points

inter_xi = -1; inter_delta = 0;

while n < len1
    [pt_xi,pt_delta, ~] = getPolarPoint(xis2, deltas2, 1, coordinate, true);
    % true because 2nd polar is an expansion
    tmp = abs(pt_xi - xis1(n));
    if tmp < distance
        inter_xi = pt_xi;
        inter_delta = pt_delta;
        distance = tmp;
    end
    n = n+1;
    coordinate = deltas1(n);
end

if inter_xi == -1
    disp("No intersection point found");
end

end