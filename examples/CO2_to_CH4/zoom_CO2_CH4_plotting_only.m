%% This script is designed to fastly get the CO2/CH4 regime diagram
% by loading the pre-calculated limits

load('sw_detonation/examples/CO2_to_CH4/CO2_CH4_poly_correction.mat')

figure
hold on

% plot the pre-calculated boundaries, uncomment lines 15 and 16 to get
% Mathilde's corrections plotted

plot(omegas_FPR_TNR,chis_FPR_TNR,'g', 'LineWidth',2,'MarkerSize',10, "DisplayName","FPR<->TNR") %ploting FPR<->TNR lim Yann's method
plot(omegas_FPR_TNR_M,chis_FPR_TNR_M,'--g', 'LineWidth',2,'MarkerSize',10, "DisplayName","FPR<->TNR (constant coefficient)") %ploting FPR<->TNR lim Mathilde's method
%plot(omegas_FPR_TNR_2,chis_FPR_TNR_2,'-.g', 'LineWidth',2,'MarkerSize',10, "DisplayName","FPR<->TNR (fitting of experimental data)") %ploting FPR<->TNR lim Mathilde's method
plot(omegas_TNR_LSR(2:end),chis_TNR_LSR(2:end),'m', 'LineWidth',2,'MarkerSize',10, "DisplayName","TNR<->LSR")  %ploting TNR<->LSR lim

cols=['k';'b';'r';'g';'m'];
lgd=['','','',"FPR<->TNR article boundary points","TNR<->LSR article boundary points"];
for i=4:no_read_lines
    hold on
    no_read_points=points(2*i-1,1);
    plot(points(2*i-1,2:no_read_points+1),points(2*i,2:no_read_points+1),[cols(i) 'o'], "DisplayName", lgd(i),'LineWidth',2.5);
    %ploting article data
end

% Data coming from Stany Gallier's CFD simulations
fpr_ang = [45.4 41.2 53.7 48.7];
fpr_khi = [0.77 0.5 0.77 0.5];
plot(fpr_ang,fpr_khi,'*g','LineWidth',2,'MarkerSize',10,"DisplayName", "FPR CFD sims");

tnr_ang = [60.8 63 50.7 52.9 60.1];
tnr_khi = [0.77 0.77 0.5 0.5 0.5];
plot(tnr_ang,tnr_khi,'*m','LineWidth',2,'MarkerSize',10,"DisplayName", "TNR CFD sims");

lsr_ang = [68.1];
lsr_khi = [0.5];
plot(lsr_ang,lsr_khi,'*k','LineWidth',2,'MarkerSize',10, "DisplayName", "LSR CFD sims");

% other settings
legend("Location","SouthEast");
set(gca,'FontSize',16,'FontWeight','bold','LineWidth',1.5)
xlabel("\omega_i (deg)",'FontSize',20, 'FontWeight', 'bold')
ylabel("\chi", 'FontSize', 20, 'FontWeight', 'bold')
xlim([48,90])
ylim([0.3,1])
hold off
%saveas(gcf,('CO2_CH4_limits_shockwaves.eps'),'epsc2');