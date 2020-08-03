function [omega_r] = solveShockAngle(M, gamma, delta)
syms om
eqn = postShockDeflection(M, gamma, om) == delta;
s = vpasolve(eqn, om);
omega_r = s;
end