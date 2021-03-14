% This script is designed to plot the dependance between alpha and omega_i
% in Stany's CFD simulations for a planar interface and a planar choc for a
% CO2/CH4 interface


figure()
%for mach=1.12 / chi = 0.77

alpha_1 = [27 32.06 33.27 34.49 38 49 65 80 88];
omega_1 = [27 32.06 33.27 34.49 38 45.4 53.7 60.8 63.0];

pente_1 = (omega_1(end)-omega_1(1))/(alpha_1(end)-alpha_1(1));
lgd_1 = strcat("chi = 0.77 and \omega = ",num2str(pente_1), "*\alpha");
plot(alpha_1,omega_1,'o-k','LineWidth',2,'MarkerSize',10, 'DisplayName', lgd_1);
hold on;

%for chi = 0.5
alpha_2 = [33 35 40 42 52 55 60 80 89.9];
omega_2 = [33 35 40 41.2 48.7 50.7 52.9 60.1 68.1];

pente_2 = (omega_2(end)-omega_2(1))/(alpha_2(end)-alpha_2(1));
lgd_2 = strcat("chi = 0.5 and \omega = ",num2str(pente_2), "*\alpha");
plot(alpha_2,omega_2,'o-b','LineWidth',2,'MarkerSize',10, 'DisplayName', lgd_2);

legend('Location', 'SouthEast');
set(gca,'FontSize',16,'FontWeight','bold','LineWidth',1.5)
xlabel("\alpha (deg)",'FontSize',20, 'FontWeight', 'bold')
ylabel("\omega_i (deg)", 'FontSize', 20, 'FontWeight', 'bold')
xlim([10,90])
ylim([10,90])