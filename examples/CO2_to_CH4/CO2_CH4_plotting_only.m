%% This script is designed to fastly get the CO2/CH4 regime diagram
% by loading the pre-calculated limits

load('sw_detonation/examples/CO2_to_CH4/CO2_CH4_limits_shockwaves.mat')

figure
hold on

% plot the pre-calculated boundaries, uncomment lines 15 and 16 to get
% Mathilde's corrections plotted

plot(omegas_RRE_RRR,chis_RRE_RRR,'k--','LineWidth',2,'MarkerSize',10) %ploting RRE<->RRR limit
plot(omegas_RRE(2:end-1),chis_RRE(2:end-1),'k', 'LineWidth',2,'MarkerSize',10) %ploting RRE<->... limit
plot(omegas_RRR_BPR,chis_RRR_BPR,'b','LineWidth',2,'MarkerSize',10) %ploting RRR<->BPR limit Yann's method
plot(omegas_BPR_FNR,chis_BPR_FNR,'r','LineWidth',2,'MarkerSize',10) %ploting BPR<->FNR lim
plot(omegas_FPR_TNR,chis_FPR_TNR,'g', 'LineWidth',2,'MarkerSize',10) %ploting FPR<->TNR lim Yann's method
plot(omegas_FPR_TNR_2,chis_FPR_TNR_2,'-.g', 'LineWidth',2,'MarkerSize',10) %ploting FPR<->TNR lim Mathilde's method
plot(omegas_TNR_LSR(2:end),chis_TNR_LSR(2:end),'m', 'LineWidth',2,'MarkerSize',10)  %ploting TNR<->LSR lim

cols=['k';'b';'r';'g';'m'];
for i=1:no_read_lines
    hold on
    no_read_points=points(2*i-1,1);
    plot(points(2*i-1,2:no_read_points+1),points(2*i,2:no_read_points+1),[cols(i) 'o'],'LineWidth',2.5);
    %ploting article data
end

vw_art_ang=[27,32.06,33.27,34.49,38,46];
n_vw_art_ang=length(vw_art_ang);
plot(vw_art_ang,.78*ones(1,n_vw_art_ang),'<k','LineWidth',2,'MarkerSize',10)
w_art_ang=[50.5,55];
n_w_art_ang=length(w_art_ang);
plot(w_art_ang,.53*ones(1,n_w_art_ang),'^k','LineWidth',2,'MarkerSize',10)
s_art_ang=[30,46,58];
n_s_art_ang=length(s_art_ang);
plot(s_art_ang,.18*ones(1,n_s_art_ang),'>k','LineWidth',2,'MarkerSize',10)

% Data coming from Stany Gallier's CFD simulations
rre_ang = [27 32.06 33];
rre_khi = [0.77 0.77 0.5];
plot(rre_ang,rre_khi,'*k','LineWidth',2,'MarkerSize',10);

rrr_ang = [33.27 34.49];
rrr_khi = [0.77 0.77];
plot(rrr_ang,rrr_khi,'*b','LineWidth',2,'MarkerSize',10);

bpr_ang = [38 35 40];
bpr_khi = [0.77 0.5 0.5];
plot(bpr_ang,bpr_khi,'*r','LineWidth',2,'MarkerSize',10);

fpr_ang = [45.4 41.2 53.7 48.7];
fpr_khi = [0.77 0.5 0.77 0.5];
plot(fpr_ang,fpr_khi,'*g','LineWidth',2,'MarkerSize',10);

tnr_ang = [60.8 63 50.7 52.9 60.1];
tnr_khi = [0.77 0.77 0.5 0.5 0.5];
plot(tnr_ang,tnr_khi,'*m','LineWidth',2,'MarkerSize',10);

lsr_ang = [68.1];
lsr_khi = [0.5];
plot(lsr_ang,lsr_khi,'ok','LineWidth',2,'MarkerSize',10);

ltxt=sprintf(...
"RRE->Intromission->RRR->Shock critical\n->BPR->FPR from Nourgaliev et al. fem sims");

% Legend without correction
% legends=legend("RRE<->RRR (analytical resolution)",...
%     "RRE->... (graphical resolution)",...
%     "RRR<->BPR","BPR<->FNR","FPR<->TNR",...
%     "TNR<->LSR",...
%     "RRE<->... article boundary points",...
%     "RRR<->BPR article boundary points",...
%     "BPR<->FNR article boundary points",...
%     "FPR<->TNR article boundary points",...
%     "TNR<->LSR article boundary points",ltxt,...
%     "TRR->TNR from Nourgaliev et al. sims",...
%     "RRE->BPR->TMR from Nourgaliev et al. fems sims",...
%     "RRE CFD sims", "RRR CFD sims","BPR CFD sims",...
%     "FPR CFD sims","TNR CFD sims","LSR CFD sims");

% Legend with correction 
legends=legend("RRE<->RRR (analytical resolution)",...
    "RRE->... (graphical resolution)",...
    "RRR<->BPR","BPR<->FNR","FPR<->TNR",...
    "FPR<->TNR (correction)","TNR<->LSR",...
    "RRE<->... article boundary points",...
    "RRR<->BPR article boundary points",...
    "BPR<->FNR article boundary points",...
    "FPR<->TNR article boundary points",...
    "TNR<->LSR article boundary points",ltxt,...
    "TRR->TNR from Nourgaliev et al. sims",...
    "RRE->BPR->TMR from Nourgaliev et al. fems sims",...
    "RRE CFD sims", "RRR CFD sims","BPR CFD sims",...
    "FPR CFD sims","TNR CFD sims","LSR CFD sims");

% other settings
legends.Location = 'NorthEast';
set(gca,'FontSize',16,'FontWeight','bold','LineWidth',1.5)
xlabel("\omega_i (deg)",'FontSize',20, 'FontWeight', 'bold')
ylabel("\chi", 'FontSize', 20, 'FontWeight', 'bold')
xlim([25,90])
ylim([0,1])
hold off
%saveas(gcf,('CO2_CH4_limits_shockwaves.eps'),'epsc2');