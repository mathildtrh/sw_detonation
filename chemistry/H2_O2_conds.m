%gas parameters
    % Molar fractions
    XH2=2/3;
    XO2=1/3;
    XHE = 1;
    % Molecular masses
    MF = 2; % fuel = H2
    MO = 32; % oxydant = O2
    MI = 4; % inert gas = He
    %names
    name_I='2/3*H2+1/3*O2';
    name_II='He';
    %ratios of specific heat
    gamma_I=1.4016; %slow material, first phase
    gamma_II=1.667; %fast material, second phase
    %molecular masses:
    mu_I=XH2*MF+XO2*MO; 
    mu_II=XHE*MI;
    %Universal gas constant (kg.mol^-1.K^-1)
    R_u=8.314;
    %Gas constants
    R_I=R_u/(mu_I*10^-3);
    R_II=R_u/(mu_II*10^-3);