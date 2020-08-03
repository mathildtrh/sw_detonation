M1 = linspace(1,10,10);
phi = linspace(0, pi);
gamma = 1.4016;
press = zeros(1,100);

figure()
hold on
for m = 1:1:10
    for g = 1:1:100
        press(g) = machtoPressJump(M1(m), gamma_I, phi(g));
    end
    plot(phi, press);
end
hold off