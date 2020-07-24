%% Examples of polar diagrams --> better understanding of transition between RRE, RRR and BPR

run H2_O2_conds.m
T1 = 600;
T5 = 1138;
a1 = sqrt(gamma_I*R_I*T1);
a5 = sqrt(gamma_II*R_II*T5);

chi = 0.3;
xi_i = 1/chi;
Msh = sqrt(xiToSqMach(xi_i, gamma_I, pi/2));

omega_deg = 26;
omega_i = omega_deg*pi/180;

M1 = Msh/sin(omega_i);
M5 = a1/a5*M1;
M2n = machtoPSNormalMach(M1, gamma_I, omega_i);
delta_1 = postShockDeflection(M1, gamma_I, omega_i);
M2 = M2n / sin(omega_i - delta_1); 
ny= 500;

figure()
[inc_x, inc_d]=getPolarMathilde(M1, gamma_I, ny, 0);
[exp_x, exp_d] = getExpPolar(M2, gamma_I, ny, xi_i, delta_1);
[refl_x, refl_d] = getPolarMathilde(M2, gamma_I, ny, 2, xi_i, delta_1);       % Polar of reflected expansion
[trans_x, trans_d] = getPolarMathilde(M5, gamma_II, ny);                       % Polar of transmitted shock

inc = semilogy(inc_d*180/pi,inc_x,'r', 'DisplayName', 'Incident shock');
hold on;
trans = semilogy(trans_d*180/pi,trans_x,'g', 'DisplayName', 'Transmitted shock');
refl = semilogy(refl_d*180/pi,refl_x,'y', 'DisplayName', 'Reflected shock');
exp = semilogy(exp_d*180/pi, exp_x, 'c', 'DisplayName', 'Reflected expansion');
%inter = semilogy(delta_t*180/pi, xi_t, '+k', 'DisplayName', 'Intersection point : (\delta_t, \xi_t)');

inc.LineWidth = 2;
trans.LineWidth = 2;
refl.LineWidth = 2;
exp.LineWidth = 2;
%inter.LineWidth = 2;

grid on
grid minor
xlabel('\omega (in deg)')
ylabel('\xi')
legend('Location', 'eastoutside')
title({'Polar figure for RRE'; strcat('\chi = ', num2str(chi),'; M1 = ', num2str(M1)); strcat('M5 = ', num2str(M5),'; \omega_i = ', num2str(omega_deg))})
