
gamma_I = 1.4016;
M2a = 4.89;
P2a = 4.5613e+05;

delta_2 = linspace(0.04,0.05);
P3 = zeros(1, length(delta_2));

for i = 1:1:length(delta_2)
    delta_3 = -0.1163 - delta_2(i);
    omega_r = solveShockAngle(M2a, gamma_I, delta_2(i));
    omega_r = double(omega_r); % if solution exists
    
    M2bn = machtoPSNormalMach(M2a, gamma_I, omega_r);
    M2b = M2bn/abs(sin(omega_r - delta_2(i)));
    P2b = P2a*machtoPressJump(M2a, gamma_I, omega_r);
    
    [~,nu_M2b, ~] = flowprandtlmeyer(gamma_I, M2b, 'mach');
    nu_M3 = delta_3 + nu_M2b;
    [M3, ~, ~] = flowprandtlmeyer(gamma_I, nu_M3, 'nu');
    
    P3(i) = P2b*pressureDropExp(M2b, M3, gamma_I);
end

plot(delta_2, P3);