%% This script is designed to fastly get the Ar/He regime diagram
% by loading the pre-calculated limits

load('/home/mathilde/Documents/PRe_Detonation/sw_detonation/plots_src/Ar_to_He/Ar_He_limits_shockwaves.mat')

% plot the pre-calculated boundaries

figure
hold on
plot(omegas_RRE(2:end-1),chis_RRE(2:end-1), 'k','LineWidth',2,'MarkerSize',10) %ploting RRE<->... limit
plot(omegas_RRR_BPR,chis_RRR_BPR,'b','LineWidth',2,'MarkerSize',10) %ploting RRR<->BPR limit
plot(omegas_BPR_NFR,chis_BPR_NFR,'r','LineWidth',2,'MarkerSize',10) %ploting BPR<->NFR lim
plot(omegas_FPR_TNR(2:end),chis_FPR_TNR(2:end),'g','LineWidth',2,'MarkerSize',10) %ploting FPR<->TNR lim
plot(omegas_TNR_LSR(2:end),chis_TNR_LSR(2:end),'m','LineWidth',2,'MarkerSize',10)  %ploting TNR<->LSR lim

Ar_to_He_sim_results_data;
plot(RRE_sims(:,3),RRE_sims(:,1),'*k','LineWidth',2,'MarkerSize',10)
plot(BPR_sims(:,3),BPR_sims(:,1),'*r','LineWidth',2,'MarkerSize',10)
plot(FPR_sims(:,3),FPR_sims(:,1),'*g','LineWidth',2,'MarkerSize',10)
plot(TMR_sims(:,3),TMR_sims(:,1),'*m','LineWidth',2,'MarkerSize',10)
legends={"RRE->... (graphical resolution)",...
    "RRR<->BPR","BPR<->FNR","FPR<->TNR",...
    "TNR<->LSR","RRE sims","BPR sims",...
    "FPR sims","TNR sims"};

% other settings
legend(legends,'Location','NorthEast')
set(gca,'FontSize',16,'FontWeight','bold','LineWidth',1.5)
xlabel("\omega_i (deg)",'FontSize',20, 'FontWeight', 'bold')
ylabel("\chi", 'FontSize', 20, 'FontWeight', 'bold')
xlim([15,90])
ylim([0,1])
hold off
saveas(gcf,('Ar_He_limits_shockwaves.eps'),'epsc2');