M2a = 4.89;
gamma_I = 1.4016;
omega_r = linspace(0,pi);
M2b = zeros(1, length(omega_r));

for i = 1:1:length(omega_r)
    M2bn = machtoPSNormalMach(M2a, gamma_I, omega_r(i));
    delta_2 = postShockDeflection(M2a, gamma_I, omega_r(i));
    M2b(i) = M2bn/sin(omega_r(i) - delta_2);
end

plot(omega_r, M2b, 'b');
hold on
plot([0 pi], [1 1], 'r');