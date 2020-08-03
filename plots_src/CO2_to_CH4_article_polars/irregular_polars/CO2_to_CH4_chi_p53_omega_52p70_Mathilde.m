%% Slight difference with article polars


%gas parameters
    % names
        name_I='CO2';
        name_II='CH4';
    % ratios of constant heat
        gamma_I=1.288;
        gamma_II=1.301;
    % molecular masses
        mu_I=44.01;
        mu_II=18.84;
    % temperature in K
        T_I = 298; 
        T_II = 298;
    % Perfect gases constant
        Ru = 8.314; 
    % speed of sound
        a_I = sqrt(gamma_I*Ru/mu_I*T_I); 
        a_II = sqrt(gamma_II*Ru/mu_II*T_II);
        
%Plotting parameters
ny=200;
%Params of figure 6a of slow-fast gas interface
chi=.53;
omega_deg=52.70;
xi_i=1/chi;
omega_rad=pi/180*omega_deg;

Msi=sqrt(xiToSqMach(xi_i,gamma_I,pi/2)); %i shock Mach
M1i=Msi/sin(omega_rad); %i free-stream Mach

Vi = Msi*a_I;

%Computing Mst:
b=(gamma_II+1)/(gamma_I+1)*(Vi^2-a_I^2)/Vi;
Vt=1/2*(b+sqrt(b^2+4*a_II^2)); 
Mst = Vt/a_II;%t shock Mach

% Let's take the membrane inertia into account
b_corrected = (gamma_II+1)/(gamma_I+1)*(Vi^2-a_I^2)/Vi * (x1*omega_rad^2 + x2*omega_rad +x3);
Vt_corrected = 1/2*(b_corrected+sqrt(b_corrected^2+4*a_II^2));
Mst_corrected = Vt_corrected/a_II;

xi_t=((1-gamma_II)+2*gamma_II*Mst^2)/(1+gamma_II); %transmitted shock...
    %... pressure jump
xi_t_corrected = ((1-gamma_II)+2*gamma_II*Mst_corrected^2)/(1+gamma_II);

%computing Msj, j shock Mach:
syms Msj
eqn= ((1-gamma_I)+2*gamma_I*Msj^2)/(1+gamma_I)*...
    (1-(gamma_II-1)/(gamma_I+1)*sqrt((gamma_I*mu_II)/(gamma_II*mu_I))*...
    (Msj-1/Msj))^(-2*gamma_II/(gamma_II-1))==xi_t;
a=[vpasolve(eqn,Msj,[1,Inf])];

syms Msj_corrected
eqn2= ((1-gamma_I)+2*gamma_I*Msj_corrected^2)/(1+gamma_I)*...
    (1-(gamma_II-1)/(gamma_I+1)*sqrt((gamma_I*mu_II)/(gamma_II*mu_I))*...
    (Msj_corrected-1/Msj_corrected))^(-2*gamma_II/(gamma_II-1))==xi_t_corrected;
a2=[vpasolve(eqn2,Msj_corrected,[1,Inf])];

if ~(isempty(a)) %checking if solution was found
    figure
    hold on
    plotPolarMathilde(M1i,gamma_I,ny); %plotting i polar

    plotPolarMathilde(Mst,gamma_II,ny); %plotting t polar
    plotPolarMathilde(Mst_corrected,gamma_II,ny); %plotting t polar
    %plotting reflected polar:
    
    Mr=sqrt(postShockMachSq(xi_i,M1i,gamma_I)); %reflected shock free-...
        %... stream Mach
    phi_i=asin(Msi/M1i);
    i_dev=postShockDeflection(M1i, gamma_I, phi_i); %i deviation
    plot([-1.1,1.1]*i_dev*180/pi,xi_i*[1,1]);
    plotPolarMathilde(Mr,gamma_I,ny,2,xi_i,i_dev); %plotting r polar

    Msj=double(a(1)); %j shock Mach
    Msj_corrected=double(a2(1));
    
    xi_j=((1-gamma_I)+2*gamma_I*Msj^2)/(1+gamma_I); %j pressure jump
    xi_j_corrected=((1-gamma_I)+2*gamma_I*Msj_corrected^2)/(1+gamma_I);
    
    M1k=sqrt(postShockMachSq(xi_j,M1i,gamma_I)); %k shock Mach
    M1k_corrected=sqrt(postShockMachSq(xi_j_corrected,M1i,gamma_I));
    
    phi_j=pi-asin(Msj/M1i); %because j shock makes an obtuse angle with incident flow
    j_dev=postShockDeflection(M1i, gamma_I, phi_j); %j deviation
    
    phi_j_corrected=pi-asin(Msj_corrected/M1i); %because j shock makes an obtuse angle with incident flow
    j_dev_corrected=postShockDeflection(M1i, gamma_I, phi_j_corrected); %j deviation
    
    plotPolarMathilde(M1k,gamma_II,ny,2,xi_j,j_dev); %plotting k polar
    plotPolarMathilde(M1k_corrected,gamma_II,ny,2,xi_j_corrected,j_dev_corrected); %plotting k polar

    
    fid=fopen('prec_pol2.txt','r');
    extract_plots;
    semilogy(points(1,2:points(1,1)+1),points(2,2:points(1,1)+1),'mo')
    semilogy(points(3,2:points(3,1)+1),points(4,2:points(3,1)+1),'go')
    
    legend("Incident and j precursor shocks","Transmitted Shock","Corrected Transmitted Shock",...
        "Incident \xi","Reflected shock", "k shock","Corrected k shock","Article r polar",...
        "Article k polar");
    title(['Polars for t-precursor wave ' name_I '->'...
         name_II ' refraction with \chi=' num2str(chi) ' and \omega='...
        num2str(omega_deg) ' deg'])
    xlabel('\delta (rad)')
    set(gca,'yscale','log')
    ylabel('\xi')
    hold off
else
    'Could not compute j shock wave Mach!'
end
