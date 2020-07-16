function [c, rpdr] = solveGeomExp(r, mu, dtheta)

% This function solves the system of equations that rules the geometrical
% evolution of the lagrangian particle when it crosses an expansion fan

plus = mu + dtheta;
minus = mu - dtheta;

syms u v
eqns = [ u == u*r*sin(plus)/(sin(minus)*(r+u*cos(plus))) - r/(cos(minus)) , v^2 == r* ( (r+u*cos(plus))*sin(plus)/sin(minus) - u*(sin(plus))^2/cos(minus))];
[c, rpdr] = vpasolve(eqns, [u v]);

c = abs(double(c(1)));
rpdr = abs(double(rpdr(1)));

end