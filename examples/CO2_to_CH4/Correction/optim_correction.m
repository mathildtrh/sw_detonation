% This MATLAB script is designed to minimize a function --> TO BE CONTINUED

% Gas parameters
gamma_I = 1.288; % specific heat ratio
gamma_II = 1.303;
mu_I = 44.01*1e-3; % molecular mass in kg/mol
mu_II = 16.04*1e-3;
T_I = 298; % temperature in K
T_II = 298;
Ru = 8.314; % Perfect gases constant
a_I = sqrt(gamma_I*Ru/mu_I*T_I); %speed of sound
a_II = sqrt(gamma_II*Ru/mu_II*T_II);

%% Choose mode
% mode 1 : polynomial
% mode 2 : cos(omega_i)
% mode 3 : integrated in calculation of boudaries

mode = 1;
if mode == 1
    %% Figure 8a to be fitted
    Ma = csvread('fig8a_hend78.csv'); % digitized fig8a
    omegas_a = Ma(:,1)*pi/180;
    VtVi_a = Ma(:,2);
    
    xi_a = 1/0.78; % Pressure jump for fig8c
    Vi_a = sqrt(xiToSqMach(xi_a, gamma_I, pi/2))*a_I; % Associated shock velocity for fig8a
    ba = (gamma_II+1)/(gamma_I+1) * (Vi_a^2 - a_I^2)/Vi_a; % Auxiliary value
    
    % Here x(1) = a; x(2) = b and x(3) = c
    % where Vpt = Vpi*(a*omega^2 + b*omega + c)
    func_a = @(x)0.5/Vi_a*(ba*(x(1)*omegas_a.^2+x(2)*omegas_a+x(3)) + sqrt((ba*(x(1)*omegas_a.^2+x(2)*omegas_a+x(3))).^2 + 4*a_II^2)) - VtVi_a;
    x0_a = [1/3 1/3 1/3];
    
    x_a = lsqnonlin(func_a, x0_a);
    
    
    %% Figure 8b to be fitted
    Mb = csvread('fig8b_hend78.csv'); % digitized fig8b
    omegas_b = Mb(:,1)*pi/180;
    VtVi_b = Mb(:,2);
    
    xi_b = 1/0.53; % Pressure jump for fig8c
    Vi_b = sqrt(xiToSqMach(xi_b, gamma_I, pi/2))*a_I; % Associated shock velocity for fig8a
    bb = (gamma_II+1)/(gamma_I+1) * (Vi_b^2 - a_I^2)/Vi_b; % Auxiliary value
    
    % Here x(1) = a; x(2) = b and x(3) = c
    % where Vpt = Vpi*(a*omega^2 + b*omega + c)
    func_b = @(x)0.5/Vi_b*(bb*(x(1)*omegas_b.^2+x(2)*omegas_b+x(3)) + sqrt((bb*(x(1)*omegas_b.^2+x(2)*omegas_b+x(3))).^2 + 4*a_II^2)) - VtVi_b;
    x0_b = [1/3 1/3 1/3];
    
    x_b = lsqnonlin(func_b, x0_b);
    
    
    % Let's start with fig8c because it is the worse
    Mc = csvread('fig8c_hend78.csv'); % digitized fig8c
    omegas_c = Mc(:,1)*pi/180;
    VtVi_c = Mc(:,2);
    
    xi_c = 1/0.18; % Pressure jump for fig8c
    Vi_c = sqrt(xiToSqMach(xi_c, gamma_I, pi/2))*a_I; % Associated shock velocity for fig8a
    bc = (gamma_II+1)/(gamma_I+1) * (Vi_c^2 - a_I^2)/Vi_c; % Auxiliary value
    
    % Here x(1) = a; x(2) = b and x(3) = c
    % where Vpt = Vpi*(a*omega^2 + b*omega + c)
    func_c = @(x)0.5/Vi_c*(bc*(x(1)*omegas_c.^2+x(2)*omegas_c+x(3)) + sqrt((bc*(x(1)*omegas_c.^2+x(2)*omegas_c+x(3))).^2 + 4*a_II^2)) - VtVi_c;
    x0_c = [1/3 1/3 1/3];
    
    x_c = lsqnonlin(func_c, x0_c);
end

if mode == 2
    
    %% Figure 8a to be fitted
    Ma = csvread('fig8a_hend78.csv'); % digitized fig8a
    omegas_a = Ma(:,1)*pi/180;
    VtVi_a = Ma(:,2);
    
    xi_a = 1/0.78; % Pressure jump for fig8c
    Vi_a = sqrt(xiToSqMach(xi_a, gamma_I, pi/2))*a_I; % Associated shock velocity for fig8a
    ba = (gamma_II+1)/(gamma_I+1) * (Vi_a^2 - a_I^2)/Vi_a; % Auxiliary value
    
    % Here x(1) = a
    % where Vpt = Vpi*a*cos(omega)
    func_a = @(x)0.5/Vi_a*(ba*x(1)*cos(omegas_a) + sqrt((ba*x(1)*cos(omegas_a)).^2 + 4*a_II^2)) - VtVi_a;
    x0_a = [1];
    
    x_a = lsqnonlin(func_a, x0_a);
    
    
    %% Figure 8b to be fitted
    Mb = csvread('fig8b_hend78.csv'); % digitized fig8b
    omegas_b = Mb(:,1)*pi/180;
    VtVi_b = Mb(:,2);
    
    xi_b = 1/0.53; % Pressure jump for fig8c
    Vi_b = sqrt(xiToSqMach(xi_b, gamma_I, pi/2))*a_I; % Associated shock velocity for fig8a
    bb = (gamma_II+1)/(gamma_I+1) * (Vi_b^2 - a_I^2)/Vi_b; % Auxiliary value
    
    % Here x(1) = a
    % where Vpt = Vpi*a*cos(omega)
    func_b = @(x)0.5/Vi_b*(bb*x(1)*cos(omegas_b) + sqrt((bb*x(1)*cos(omegas_b)).^2 + 4*a_II^2)) - VtVi_b;
    x0_b = [1];
    
    x_b = lsqnonlin(func_b, x0_b);
    
    
    Mc = csvread('fig8c_hend78.csv'); % digitized fig8c
    omegas_c = Mc(:,1)*pi/180;
    VtVi_c = Mc(:,2);
    
    xi_c = 1/0.18; % Pressure jump for fig8c
    Vi_c = sqrt(xiToSqMach(xi_c, gamma_I, pi/2))*a_I; % Associated shock velocity for fig8a
    bc = (gamma_II+1)/(gamma_I+1) * (Vi_c^2 - a_I^2)/Vi_c; % Auxiliary value
    
    % Here x(1) = a
    % where Vpt = Vpi*a*cos(omega)
    func_c = @(x)0.5/Vi_c*(bc*x(1)*cos(omegas_c) + sqrt((bc*x(1)*cos(omegas_c)).^2 + 4*a_II^2)) - VtVi_c;
    x0_c = [1];
    
    x_c = lsqnonlin(func_c, x0_c);
end

if mode == 3
    %% Figure 8b to be fitted
    Mb = csvread('fig8b_hend78.csv'); % digitized fig8b
    omegas_b = Mb(:,1)*pi/180;
    VtVi_b = Mb(:,2);
    
    xi_b = 1/0.53; % Pressure jump for fig8c
    Vi_b = sqrt(xiToSqMach(xi_b, gamma_I, pi/2))*a_I; % Associated shock velocity for fig8a
    bb = (gamma_II+1)/(gamma_I+1) * (Vi_b^2 - a_I^2)/Vi_b; % Auxiliary value
    
    % Here x(1) = a; x(2) = b and x(3) = c
    % where Vpt = Vpi*(a*cos(omega) + b*sin(omega) + c)
    func_b = @(x)0.5/Vi_b*(bb*(x(1)*cos(x(4)*omegas_b)+x(2)*sin(x(5)*omegas_b)+x(3)) + sqrt((bb*(x(1)*cos(x(4)*omegas_b)+x(2)*sin(x(5)*omegas_b)+x(3))).^2 + 4*a_II^2)) - VtVi_b;
    x0_b = [1/3 1/3 1/3 1 1];
    
    x_b = lsqnonlin(func_b, x0_b);
end

if mode == 4
    %% Figure 8b to be fitted
    Mb = csvread('fig8b_hend78.csv'); % digitized fig8b
    omegas_b = Mb(:,1)*pi/180;
    VtVi_b = Mb(:,2);
    
    xi_b = 1/0.53; % Pressure jump for fig8c
    Vi_b = sqrt(xiToSqMach(xi_b, gamma_I, pi/2))*a_I; % Associated shock velocity for fig8a
    bb = (gamma_II+1)/(gamma_I+1) * (Vi_b^2 - a_I^2)/Vi_b; % Auxiliary value
    
    % Here x(1) = a; x(2) = b and x(3) = c
    % where Vpt = Vpi*(a*cos(omega) + b*sin(omega) + c)
    func_b = @(x)0.5/Vi_b*(bb*(x(1)*omegas_b.^2+x(2)*omegas_b+x(3)) + sqrt((bb*(x(1)*omegas_b.^2+x(2)*omegas_b+x(3))).^2 + 4*a_II^2)) - VtVi_b;
    x0_b = [1/3 1/3 1/3];
    
    x_b = lsqnonlin(func_b, x0_b);
end



