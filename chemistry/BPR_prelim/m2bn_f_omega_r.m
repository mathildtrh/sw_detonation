M2a = 4.89;
gamma_I = 1.4016;
omega_r = linspace(0,pi);
M2bn = zeros(1, length(omega_r));

for i = 1:1:length(omega_r)
    M2bn(i) = machtoPSNormalMach(M2a, gamma_I, omega_r(i));
end

plot(omega_r, M2bn, 'b');
hold on
plot([0 pi], [1 1], 'r');