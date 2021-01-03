% This MATLAB script is designed to get a better adjustment of FPR to TNR
% boundary by fitting Abd el Fattah and Henderson 1978 experimental
% results and to compare interpolation with their piston theory

plotting_mode = 2;

if plotting_mode == 1
    gamma_I = 1.288; % specific heat ratio
    gamma_II = 1.303;
    mu_I = 44.01*1e-3; % molecular mass in kg/mol
    mu_II = 16.04*1e-3;
    T_I = 298; % temperature in K
    T_II = 298;
    Ru = 8.314; % Perfect gases constant
    a_I = sqrt(gamma_I*Ru/mu_I*T_I); %speed of sound
    a_II = sqrt(gamma_II*Ru/mu_II*T_II);
    
    
    Ma = csvread('fig8a_hend78.csv'); % digitized fig8a
    Mb = csvread('fig8b_hend78.csv'); % digitized fig8b
    Mc = csvread('fig8c_hend78.csv'); % digitized fig8c
    
    xi_a = 1/0.78; % Pressure jump for fig8a
    xi_b = 1/0.53; % Pressure jump for fig8b
    xi_c = 1/0.18; % Pressure jump for fig8c
    
    Vi_a = sqrt(xiToSqMach(xi_a, gamma_I, pi/2))*a_I; % Associated shock velocity for fig8a
    Vi_b = sqrt(xiToSqMach(xi_b, gamma_I, pi/2))*a_I; % Associated shock velocity for fig8a
    Vi_c = sqrt(xiToSqMach(xi_c, gamma_I, pi/2))*a_I; % Associated shock velocity for fig8a
    
    ba = (gamma_II+1)/(gamma_I+1) * (Vi_a^2 - a_I^2)/Vi_a;
    bb = (gamma_II+1)/(gamma_I+1) * (Vi_b^2 - a_I^2)/Vi_b;
    bc = (gamma_II+1)/(gamma_I+1) * (Vi_c^2 - a_I^2)/Vi_c;
    
    Vt_a_th = 0.5*(ba + sqrt(ba^2 + 4*a_II^2));
    Vt_b_th = 0.5*(bb + sqrt(bb^2 + 4*a_II^2));
    Vt_c_th = 0.5*(bc + sqrt(bc^2 + 4*a_II^2));
    
    Vt_Vi_a = 1.477; % Theoretical ratio digitized from fig8a
    Vt_Vi_b = 1.397; % Theoretical ratio digitized from fig8b
    Vt_Vi_c = 1.231; % Theoretical ratio digitized from fig8c
    
    %% Plotting data
    figure(1)
    hold on
    plot(Ma(:,1),Ma(:,2)*Vi_a, 'bo', 'LineWidth', 1.5)
    %plot([20 90], [1 1]*Vt_a_th, 'y--', 'LineWidth', 1.5)
    plot([20 90], [1 1]*Vt_Vi_a*Vi_a, 'c--', 'LineWidth', 1.5)
    ylim([Vi_a 1.6*Vi_a])
    
    figure(2)
    hold on
    plot(Mb(:,1),Mb(:,2)*Vi_b, 'bo', 'LineWidth', 1.5)
    %plot([20 90], [1 1]*Vt_b_th, 'y--', 'LineWidth', 1.5)
    plot([20 90], [1 1]*Vt_Vi_b*Vi_b, 'c--', 'LineWidth', 1.5)
    ylim([Vi_b 1.6*Vi_b])
    
    figure(3)
    hold on
    plot(Mc(:,1),Mc(:,2)*Vi_c, 'bo', 'LineWidth', 1.5)
    %plot([20 90], [1 1]*Vt_c_th, 'y--', 'LineWidth', 1.5)
    plot([20 90], [1 1]*Vt_Vi_c*Vi_c, 'c--', 'LineWidth', 1.5)
    ylim([Vi_c 1.6*Vi_c])
end

if plotting_mode == 2
    
    %optimization calculus
    mode = 4;
    run optim_correction.m
    x1 = x_b(1); x2 = x_b(2); x3 = x_b(3);
    omega_deg = linspace(20, 90);
    omega_rad = omega_deg*pi/180;
    b=(gamma_II+1)/(gamma_I+1)*(Vi_b^2-a_I^2)/Vi_b *(x1*omega_rad.^2 +x2*omega_rad +x3);
    Vt_opti=1/2*(b+sqrt(b.^2+4*a_II^2));
    
    %experimental results
    Vt_Vi_b = 1.397;
    
    %plotting
    figure()
    hold on
    plot(Mb(:,1),Mb(:,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2)
    plot([20 90], [1 1]*Vt_Vi_b, 'k--', 'LineWidth', 2)
    plot(omega_deg, Vt_opti/Vi_b, 'r', 'LineWidth', 2)
    lgd=legend('Digitized experimental data',...
        'Theoretical value',...
        'Optimized ratio','Location','NorthEast');
    xlabel('\omega (deg)','FontSize',20,'FontWeight','bold')
    ylabel('V_t/V_i','FontSize',20,'FontWeight','bold')
    ylim([1 1.6])
    set(gca,'FontSize',16,'FontWeight','bold','LineWidth',1.5)
    
end



