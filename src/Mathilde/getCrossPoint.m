% This program is designed to get the intersection point between two polars
% USEFUL : getPolarPoint.m returns xi=-1 if there is no point coresponding
% to a certain set of coordinates

function [inter_delta, inter_xi] = getCrossPoint(deltast, xist, deltasr, xisr, varargin)

% mode = 1, reflected polar is an expansion
% mode = 2, reflected polar is a shock

% while no intersection point is found
% choose a point in the 1st polar (transmitted) = coordinates
% try to find a point in the 2nd polar (reflected) that reaches coordinates
% if there isn't try another point
% if the distance between the meeting point and the coordinate is minimal
% keep this point in memory

mode = varargin{1};

if mode == 1
    len1 = length(deltast);
    
    n = 1;
    coordinate = deltast(n);
    distance = max(xist); % init of distance between to points
    
    inter_xi = -1; inter_delta = 0;
    
    while n < len1
        [pt_xi,pt_delta, ~] = getPolarPoint(xisr, deltasr, 1, coordinate, true);
        % true because 2nd polar is an expansion
        tmp = abs(pt_xi - xist(n));
        if tmp < distance
            inter_xi = pt_xi;
            inter_delta = pt_delta;
            distance = tmp;
        end
        n = n+1;
        coordinate = deltast(n);
    end
    
    if inter_xi == -1
        disp("No intersection point found");
    end
    
    
else
    len1 = length(deltast);
    
    n = 1;
    coordinate = deltast(n);
    distance = max(xist); % init of distance between to points
    
    inter_xi = -1; inter_delta = 0;
    
    while n < len1
        [pt_2xi,pt_2delta, ~] = getPolarPoint(xisr, deltasr, 1, coordinate, false);
        % false because 2nd polar is a shock --> returns 2 values
        % (strong and weak shock) choose weak shock
%         if pt_2delta(1)*pt_2delta(2) < 0
%             [pt_delta,k] = max(pt_2delta);
%             pt_xi = pt_2xi(k);
%         else
%             [pt_xi, k] = min(pt_2xi);
%             pt_delta = pt_2delta(k);
%         end
        
        if pt_2delta(1)< 0
            if pt_2delta(2)< 0
                n = n+1;
                coordinate = deltast(n);
                continue
            else
            pt_delta = pt_2delta(2);
            pt_xi = pt_2xi(2);
            end
        elseif pt_2delta(2)< 0
            pt_delta = pt_2delta(1);
            pt_xi = pt_2xi(1);
        else
            [pt_delta, k] = max(pt_2delta);
            pt_xi = pt_2xi(k);
        end
        
        tmp = abs(pt_xi - xist(n));
        if tmp < distance
            inter_xi = pt_xi;
            inter_delta = pt_delta;
            distance = tmp;
        end
        n = n+1;
        coordinate = deltast(n);
    end
    
    if inter_xi == -1
        disp("No intersection point found");
    end
end
end