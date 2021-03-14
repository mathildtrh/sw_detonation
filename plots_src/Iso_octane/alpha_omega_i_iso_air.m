% This script is designed to plot the dependance between alpha and omega_i
% in Stany's CFD simulations for a curved interface and a planar choc for a
% iso-octane/air interface


figure()
%for mach=1.1 / chi = 0.8092

alpha_1 = [11.1 20.7 27.4 32.3 37.4 42.1 46.8 51.4 65.7 85.8];
omega_1 = [11.1 20.7 27.4 32.3 37.4 42.1 46.2 50.6 64.4 85.4];

pente_1 = (omega_1(end)-omega_1(1))/(alpha_1(end)-alpha_1(1));
lgd_1 = strcat("chi = 0.8092 and \omega = ",num2str(pente_1), "*\alpha");
plot(alpha_1,omega_1,'*-k','LineWidth',2,'MarkerSize',10, 'DisplayName', lgd_1);
hold on;

%for mach=1.5 / chi = 0.4161
alpha_2 = [17.4 26.9 34.3 40.4 46.2 52.6 58.1 63.9];
omega_2 = [17.4 26.9 34.3 40.4 46.2 52.0 57.2 62.8];

pente_2 = (omega_2(end)-omega_2(1))/(alpha_2(end)-alpha_2(1));
lgd_2 = strcat("chi = 0.4161 and \omega = ",num2str(pente_2), "*\alpha");
plot(alpha_2,omega_2,'*-b','LineWidth',2,'MarkerSize',10, 'DisplayName', lgd_2);

%for mach=2 / chi = 0.2290
alpha_3 = [24.3 34.4 42.5 49.6 57.1];
omega_3 = [24.3 34.4 43.5 49.6 56.4];

pente_3 = (omega_3(end)-omega_3(1))/(alpha_3(end)-alpha_3(1));
lgd_3 = strcat("chi = 0.2290 and \omega = ",num2str(pente_3), "*\alpha");
plot(alpha_3,omega_3,'*-r','LineWidth',2,'MarkerSize',10, 'DisplayName', lgd_3);

legend('Location', 'SouthEast');
set(gca,'FontSize',16,'FontWeight','bold','LineWidth',1.5)
xlabel("\alpha (deg)",'FontSize',20, 'FontWeight', 'bold')
ylabel("\omega_i", 'FontSize', 20, 'FontWeight', 'bold')
xlim([10,90])
ylim([10,90])