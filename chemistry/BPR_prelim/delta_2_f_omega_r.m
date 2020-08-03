gamma_I = 1.4016;
mach = 4.89;
omega_r = linspace(0, pi);
deltas = zeros(1, length(omega_r)-1);
for i = 1:1:length(omega_r)-1
    deltas(i) = postShockDeflection(mach, gamma_I, omega_r(i));
end

plot(omega_r(1:end-1), deltas, 'b')
hold on
plot([0 pi], [0 0], 'r')