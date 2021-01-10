%Iso_air regime diagram plotting only

load('sw_detonation/plots_src/Iso_octane/iso_air_limits_shockwaves.mat')


% Ploting limits
figure
hold on
plot(omegas_RRE(2:end-1),chis_RRE(2:end-1),'k', 'LineWidth',2,'MarkerSize',10) %ploting RRE<->... limit
plot(omegas_RRR_BPR,chis_RRR_BPR,'b','LineWidth',2,'MarkerSize',10) %ploting RRR<->BPR limit Yann's method
plot(omegas_BPR_FNR,chis_BPR_FNR,'r','LineWidth',2,'MarkerSize',10) %ploting BPR<->FNR lim
plot(omegas_FPR_TNR,chis_FPR_TNR,'g', 'LineWidth',2,'MarkerSize',10) %ploting FPR<->TNR lim Yann's method
plot(omegas_TNR_LSR(2:end),chis_TNR_LSR(2:end),'m', 'LineWidth',2,'MarkerSize',10)  %ploting TNR<->LSR lim


% Points of interest
omegas = [10 20 30 40 50 60 70 80 90];
chi_1 = [0.228 0.228 0.228 0.228 0.228 0.228 0.228 0.228 0.228];
chi_2 = [0.416 0.416 0.416 0.416 0.416 0.416 0.416 0.416 0.416];
chi_3 = [0.809 0.809 0.809 0.809 0.809 0.809 0.809 0.809 0.809];

plot(omegas, chi_1, "*k", 'LineWidth',2,'MarkerSize',10); 
plot(omegas, chi_2, "+k", 'LineWidth',2,'MarkerSize',10);
plot(omegas, chi_3, "ok", 'LineWidth',2,'MarkerSize',10);

legends=legend("RRE->... (graphical resolution)",...
    "RRR<->BPR","BPR<->FNR","FPR<->TNR","TNR<->LSR",...
    "points for \chi = 0.228", "points for \chi = 0.416", "points for \chi = 0.809");

legends.Location = 'NorthEast';
set(gca,'FontSize',16,'FontWeight','bold','LineWidth',1.5)
xlabel("\omega_i (deg)",'FontSize',20, 'FontWeight', 'bold')
ylabel("\chi", 'FontSize', 20, 'FontWeight', 'bold')
xlim([25,90])
ylim([0,1])
hold off