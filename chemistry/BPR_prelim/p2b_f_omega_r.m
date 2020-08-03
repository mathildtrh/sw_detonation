P2a = 4.56e+05;
M2a = 4.89;
gamma_I = 1.4016;
omega_r = linspace(0,pi);
P2b = zeros(1, length(omega_r));

for i = 1:1:length(omega_r)
    P2b(i) = P2a*machtoPressJump(M2a, gamma_I, omega_r(i));
end

plot(omega_r, P2b, 'b');
hold on
plot([0 pi], [P2a P2a], 'r');