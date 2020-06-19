%gas parameters
    %names
    name_I='2/3*H2+1/3*O2';
    name_II='He';
    %ratios of specific heat
    gamma_I=1.4016; %slow material, first phase
    gamma_II=1.667; %fast material, second phase
    %molecular masses:
    mu_I=2/3*2+1/3*16;
    mu_II=4;
    %Universal gas constant (kg.mol^-1.K^-1)
    R_u=8.314;
    %Gas constants
    R_I=R_u/(mu_I*10^-3);
    R_II=R_u/(mu_II*10^-3);
    
%Experimental conditions:
    %Angle of incidence:
    omega_deg=18;
    omega_rad=omega_deg*pi/180;
    %Incident wave pressure jump
        %Cond Remy
            xi=18.5;

    %inital temperature:
        %cond Remy
        T0_I=600; %Phase I
        T0_II=1138; %Phase II
        
    %initial pressure:
        %Cond Remy
        P0=101e3;

    %initial specific volume:
    Vm0_I=R_I*T0_I/P0;
    %distance of shock to interface at t=0 (m):
    experiment_dim=1e-2;
    
    %Speed of sound
    a0_I=sqrt(gamma_I*R_I*T0_I);
    a0_II=sqrt(gamma_II*R_II*T0_II);