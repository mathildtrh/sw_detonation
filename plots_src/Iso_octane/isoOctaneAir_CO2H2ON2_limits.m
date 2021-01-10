%% regime diagram for iso-octane-air // CO2-H2O-N2 system
    %%...for iso-octane-air->CO2-H2O-N2 slow-fast refraction  in...
    %%... interface angle- incident pressure jump plane.
    
limits_computation_time=cputime;%initialising timer for computation time


%% gas parameters
    %ratios of specific heat
    gamma_iso=1.2791; %slow material, first phase
    gamma_air=1.2791; %fast material, second phase
    %molecular masses:
    mu_iso=30.26195;
    mu_air=28.606992;
    %temperature ratio
    temp_ratio=820/2500;

    
%% Computing RRE<->... limit with graphical method:
    %Principle: Doing a dichotomy on "isSfStrongRRE" function to find omega...
        %at given chi.
        %"isSfStrongRRE" uses graphical method to determine if...
        %..."RRE->..." transition was made
    disp('RRE-... graphical running')
    nchis_RRE=150;nys_RRE=100;%nchis_RRE: no. points on chi-axis
        %nys_RRE: no. points on xi axis when computing polar
    chis_RRE=linspace(0,1,nchis_RRE);%chi coordinates
    omegas_RRE=zeros(1,nchis_RRE);%pre-allocated omega degree coordinates
    rre_bound_res=.01;%omega coordinate precision, dichotomy termination threshold

    %executing dichotomy for each chi value
    for i=2:nchis_RRE-1 %not solving for chis_RRE(1)=0 and chis_RRE(end)=1 limits
        omega_inf_rre=0;omega_sup_rre=90;
        while omega_sup_rre-omega_inf_rre>rre_bound_res
            next_omega=(omega_sup_rre+omega_inf_rre)/2;
            is_sf_RRE=isSfStrongRRE(1/chis_RRE(i),pi/180*next_omega,gamma_iso,...
                gamma_air,mu_iso,mu_air,nys_RRE,temp_ratio);
            if is_sf_RRE
                omega_inf_rre=next_omega;
            else
                omega_sup_rre=next_omega;
            end
        end
        omegas_RRE(i)=omega_inf_rre;
    end
  disp('RRE-... graphical done')
    
 
%% Computing RRR<->BPR limit in omega-xi plane: Yann's method
disp('RRR-BPR (Yann) running')
    nchis=100;nys=400;%nchis: number of points along chi axis;...
        %nys: nb of points along xi axis when computing polars

    %Principle: For each value of chi, the omega coordinate is computing by...
        %...doing a dichotomy with the isBeforeSfRRRToBPR, which determines if the...
        %... RRR->BPR transition has taken place at given chi and omega
    rrr_bpr_boundary_precision=.01; %omega precision, dichtomy termination...
        %...threshold

    chis_RRR_BPR=linspace(.42,1,nchis);%chi coordinates
    omegas_RRR_BPR=zeros(1,nchis);%pre-allocated omega coordinates
    for i=1:nchis
        omega_inf=0;omega_sup=90;%initial dichotomy bounds
        while omega_sup-omega_inf>rrr_bpr_boundary_precision %dichotomy...
            %...termination criteria

            %executing dichotomy, with lower bound having true value and...
                %... upper bound having false value
            next_omega=(omega_sup+omega_inf)/2;
            if isBeforeSfRRRToBPR(1/chis_RRR_BPR(i),pi/180*next_omega,gamma_iso,...
                    gamma_air,mu_iso,mu_air,nys,temp_ratio)
                omega_inf=next_omega;
            else
                omega_sup=next_omega;
            end
        end
        omegas_RRR_BPR(i)=omega_inf; %adding omega coordinate in degrees to...
            %...preallocated array
    end
disp('RRR-BPR (Yann) done')
    
   
%% Computing BPR<->FNR boundary
    %Principle: Computes omegas for which Mt=1 at different values of chi
    
    disp('BPR-FNR running')
    ny_BPR_FNR=100;%number of points on chi axis
    chis_BPR_FNR=linspace(0,1,ny_BPR_FNR);%chi coordinates
    omegas_BPR_FNR=zeros(1,ny_BPR_FNR);%pre-allocated omega coordinates
    %calculating omega coordinates
    for i=1:ny_BPR_FNR
        omegas_BPR_FNR(i)=bPRFNROmega(1/chis_BPR_FNR(i),gamma_iso,gamma_air...
            ,mu_iso,mu_air,temp_ratio);
    end
disp('BPR-FNR done')
    
    
    
%% Computes TNR<->LSR limit
    %Principle: Solves the Mr=1 equation
    disp('TNR-LSR running')
    npts=400; %number of points on chi axis
    chis_TNR_LSR=linspace(0,1,npts);
    omegas_TNR_LSR=zeros(1,npts);
    syms omega
    for i=2:npts
        %solving Mr=1 equation at specific omega
        xi=1/chis_TNR_LSR(i);
        eqn = postShockMachSq(xi,sqrt(xiToSqMach(xi,gamma_iso,omega)),gamma_iso)==1;
        a=[vpasolve(eqn,omega,[0,pi/2])];
        if ~(isempty(a))
            omegas_TNR_LSR(i)=180/pi*a;
        end
    end
disp('TNR-LSR done')
    
    
    
%% Computing FPR<->TNR boundary (Yann's method)
%Principle: Graphically resolves the beta1=beta2 problem, as described...
    %... in Abd-El-Fattah, L. F. Henderson: Shock waves at a slow-fast...
    %... gas interface (1978), p88, figure 6b. The boundary is computed...
    %... with a dichotomy on isBeforeFPRToTNR at constant chi
disp("FPR-TNR (Yann) running")
nchis_tnr=100; %number of points on chi axis
nys_tnr=30; %number of points in computed polars
chis_FPR_TNR=linspace(.1,1,nchis_tnr); %chi values
omegas_FPR_TNR=zeros(1,nchis_tnr);%pre-allocated omega array
Msj=0; %j-wave Mach initialization
fpr_tnr_boundary_precision=.01; %dichotomy termination threshold
for i=1:nchis_tnr
    solve_Msj_eq=true;
    omega_inf_tnr=0;omega_sup_tnr=90;
    while omega_sup_tnr-omega_inf_tnr>fpr_tnr_boundary_precision
        next_omega=(omega_sup_tnr+omega_inf_tnr)/2;
        [is_before_trans,Msj]=isBeforeFPRToTNR(1/chis_FPR_TNR(i),pi/180*next_omega,...
             gamma_iso,gamma_air,mu_iso,mu_air,nys_tnr,solve_Msj_eq,Msj,temp_ratio);
        if solve_Msj_eq %Msj is only computed at the beginning...
                %... of the dichotomy
            solve_Msj_eq=false;
        end
        if is_before_trans
            omega_inf_tnr=next_omega;
        else
            omega_sup_tnr=next_omega;
        end
    end
    omegas_FPR_TNR(i)=omega_inf_tnr;
end
disp('FPR-TNR (Yann) done')



limits_computation_time=cputime-limits_computation_time; %stopping timer
disp(limits_computation_time);





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
%saveas(gcf,('CO2_CH4_limits_shockwaves.eps'),'epsc2');


save('iso_air_limits_shockwaves.mat');