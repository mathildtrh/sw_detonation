% This programm computes the evolution of pressure, temperature and
% specific volume of a Lagrangian particle when it travels through a RRE
% structure

% Inputs : mode = 1 : H2_O2 // He interface
%          x, y coordinates of the lag. particle at the beginning of the
%          simulation
%          mach : mach number of the incident shock
%          plot : if true, plots the evolution of each variable
% Outputs : P1 T1 V1 pressure, temperature and spec volume after incident shock
%           P2 T2 V2 pressure, temperature and spec volume after reflected shock

function [P0, T0_I, t_0, P1, T1, t_1, P2, T2, t_2]=RRRevolution(mode, x, y, mach, plt)
if mode ~= 1
    disp("Unexisting mode")
else
    %% Initial data in zones 0I et 0II
    run H2_O2_conds.m
    r_0 = y - x/tan(omega_rad);
    if r_0 < 0
        disp("This code only follows particles in phase I!")
    else
        
        % time spent in zone 0I
        t_0 = (experiment_dim + x) / (mach*a0_I);
        %% Computation zone 1
        
        % Pressure jump
        xi = machtoPressJump(mach, gamma_I, omega_rad);
        % Pressure
        P1 = P0*xi;
        % Temperature
        T1 = T0_I*machtoTempJump(mach, gamma_I, omega_rad);
        % Volume
        V1 = R_I*T1/P1;
        % New normal mach number
        M1n = machtoPSNormalMach(mach, gamma_I, omega_rad);
        % Deflection angle
        delta_1 = postShockDeflection(mach, gamma_I, omega_rad);
        % New Mach number
        M1 = M1n / sin(omega_rad-delta_1);
        
        %% Polar computations
        ny = 500; % number of points used to calculate polars
        M0II = a0_I/a0_II*mach; % Transmitted mach number
        % Transmitted polar
        [trans_xi, trans_delta] = getPolar(M0II, gamma_II,ny,0); % Transmitted shock polar
        % Reflected polar
        [refl_xi, refl_delta] = getPolar(M1, gamma_I, ny, 2, xi, delta_1); %mode 2 = second deviation
        % Intersection point
        [delta_t, xi_t] = getCrossPoint(real(trans_delta), trans_xi,refl_delta, refl_xi);
        % If no intersection point --> wrong refraction structure
        if xi_t == -1
            "No intersection point : maybe you're not in a RRR structure..."
            return
        end
        % Plot polars if needed
        %[i_xi, i_delta]=getPolar(M0, gamma_I, ny, 0);
        %semilogy(i_delta,i_xi,'b');
        %hold on;
        %semilogy(t_delta,t_xi,'r');
        %semilogy(r_delta,r_xi,'y');
        
        %% Computation zone 2
        
        % Pressure P2 = Pt
        P2 = P0*xi_t;
        % Reflected pressure jump
        xi_r = P2/P1;
        % Normal Mach number
        M1n_prime = sqrt( (xi_r*gamma_I + xi_r + gamma_I -1)/(2*gamma_I) );
        % New normal mach number
        M2n = machtoPSNormalMach(M1n_prime, gamma_I, pi/2); %little trick with pi/2, see function
        % Temperature
        T2 = T1*machtoTempJump(M1n_prime, gamma_I, pi/2);
        % Volume
        V2 = R_I*T2/P2;
        % Deflection angle
        delta_2 = delta_t - delta_1;
        
        % angle of the reflected shock
        beta = asin(M1n_prime/M1);
        %time spent in zone 1
        t_1 = r_0*cos(omega_rad+delta_1)*(1+tan(omega_rad+delta_1)/tan(beta))/(M1*a0_I);
        %time spent in zone 2, arbitrary set to be 2ms
        t_2 = 2e-6;
        
        if plt
            % Time array
            t_plot = [0, t_0, t_0, t_0+t_1, t_0+t_1, t_0+t_1+t_2];
            
            % Pressure array
            P_plot = [P0, P0, P1, P1, P2, P2];
            
            % Temperature array
            T_plot = [T0_I, T0_I, T1, T1, T2, T2];
            
            % Spec volume array
            V_plot = [Vm0_I, Vm0_I, V1, V1, V2, V2];
            
            subplot(3,1,1)
            plot(t_plot, P_plot, '-.r')
            title('Pressure evolution')
            xlabel("Time (s)")
            ylabel("Pressure (Pa)")
            ylim([0.9e5,2.6e5])
            
            subplot(3,1,2)
            plot(t_plot, T_plot, '-.r')
            title('Temperature evolution')
            xlabel("Time (s)")
            ylabel("Temperature (K)")
            ylim([500,3000])
            
            subplot(3,1,3)
            plot(t_plot, V_plot, '-.r')
            title('Specific volume evolution')
            xlabel("Time (s)")
            ylabel("Specific volume")
            ylim([6.5,15])
        end
    end
end
end