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

function [solution]=RRRevol(cst, x, y, mach, omega_deg, press, temp_I, temp_II, varargin)

solution = -1;
gamma_I = cst(1);
gamma_II = cst(2);
R_I = cst(3);
R_II = cst(4);

% Varargin holds the plt condition and the color of the plot if needed
plt = varargin{1};

% Is the particle in the right phase (I)?
if y < 0
    disp("Wrong phase, change coordinates. EXIT")
else
    omega_i = omega_deg*pi/180; % angle of incidence in rad
    % Is the particle before the incident shock?
    if x - y/tan(omega_i) < 0
        disp("Already behind the shock, change coordinates. EXIT")
    else
         %% Zone 1
        
        P1 = press;                                                            % Pressure
        T1 = temp_I;                                                           % Temperature
        V1 = R_I*T1/P1;                                                        % Spec. volume
        a1 = sqrt(gamma_I*R_I*T1);                                             % Speed of sound
        
        %% Change of frame of reference
        
        M1 = mach/sin(omega_i);                                                % Mach number in zone 1
        u1 = M1*a1;                                                            % u1 = -u_node -> change of frame of reference
        
        t_1 = 1/u1 * (x - y/tan(omega_i));                                     % Time spent in zone 1
        
        %% Zone 5
        
        P5 = press;                                                            % Pressure
        T5 = temp_II;                                                          % Temperature
        V5 = R_II*T5/P5;                                                       % Spec. volume
        a5 = sqrt(gamma_II*R_II*T5);                                           % Speed of sound
        M5 = a1/a5*M1;                                                         % Mach number
        
        %% Zone 2 (oblique shock from 1)
        
        M2n = machtoPSNormalMach(M1, gamma_I, omega_i);                        % Mach number (2) normal to incident shock
        xi_i = machtoPressJump(M1, gamma_I, omega_i);                          % Pressure jump
        P2 = xi_i*P1;                                                          % Pressure
        V2 = V1*machtoVolJump(M1, gamma_I, omega_i);                           % Spec. volume
        T2 = T1*machtoTempJump(M1, gamma_I, omega_i);                          % Temperature
        a2 = sqrt(gamma_I*R_I*T2);                                             % Speed of sound
        delta_1 = postShockDeflection(M1, gamma_I, omega_i);                   % Deflection angle behind the shock
        
        M2 = M2n / sin(omega_i - delta_1);                                     % Mach number
        u2 = M2*a2;                                                            % Velocity of flow
        
        %% Zone 4 (oblique shock from 5)
        ny = 500;                                                              % Number of points used to calculate polars
        
        [refl_xi, refl_delta] = getPolarMathilde(M2, gamma_I, ny, 2, xi_i, delta_1);   % Reflected shock polar : mode 2 = second deviation
        [trans_xi, trans_delta] = getPolarMathilde(M5, gamma_II,ny,0);                 % Transmitted shock polar
        
%         if ~ isreal(refl_delta)
%             disp("Complex deviation : maybe you're not in a RRR structure...");
%             return
%         end
        
        [delta_t, xi_t] = getCrossPoint(trans_delta, trans_xi,refl_delta, refl_xi,2);        
        % If no intersection point --> wrong refraction structure
        
        if xi_t == -1
            disp("No intersection point : maybe you're not in a RRR structure...");
            return
        end
        
%         % Plot polars if needed
%         figure()
%         [inc_xi, inc_delta]=getPolarMathilde(M1, gamma_I, ny, 0);
%         %[exp_xi, exp_delta] = getExpPolar(M2, gamma_I, ny, xi_i, delta_1);      
%         inc = semilogy(inc_delta*180/pi,inc_xi,'r', 'DisplayName', 'Incident shock');
%         hold on;
%         %exp = semilogy(exp_delta*180/pi,exp_xi,'c', 'DisplayName','Reflected expansion');
%         trans = semilogy(trans_delta*180/pi,trans_xi,'g', 'DisplayName', 'Transmitted shock');
%         refl = semilogy(refl_delta*180/pi,refl_xi,'y', 'DisplayName','Reflected shock');
%         inter = semilogy([delta_t*180/pi],[xi_t],'+k', 'DisplayName','Intersection point : (\delta_t, \xi_t)');
%         inter.LineWidth = 2;
%         inc.LineWidth = 2;
%         %exp.LineWidth = 2;
%         trans.LineWidth = 2;
%         refl.LineWidth = 2;
%         xlabel('\omega (in deg)')
%         ylabel('\xi')
%         legend('Location', 'eastoutside')
%         grid on
%         grid minor
%         title({'Polar figure for RRR'; strcat('M1 = ', num2str(M1), '; Particle x = ', num2str(x*100))})

        P4 = xi_t*P5;                                                          % Pressure
        V4 = V5*xiToVolJump(xi_t, gamma_II);                                   % Spec. Volume
        T4 = T5*xi_t*xiToVolJump(xi_t, gamma_II);                              % Temperature
        a4 = sqrt(gamma_II*R_II*T4);                                           % Speed of sound
        omega_t = asin(sqrt(xiToSqMach(xi_t, gamma_II, pi/2))/M5);             % omega_t = asin(M5n/M5)

        %% Zone 3 (oblique shock from 2)
        
        P3 = P4;                                                               % Pressure P3 = P4 = P5*xi_t
        xi_r = P3/P2;                                                          % Reflected pressure jump
        M2n_prime = sqrt( xiToSqMach(xi_r, gamma_I, pi/2) );                   % Normal Mach number relative to reflected shock
        M3n = machtoPSNormalMach(M2n_prime, gamma_I, pi/2);                    % New normal mach number = little trick with pi/2, see function
        T3 = T2*machtoTempJump(M2n_prime, gamma_I, pi/2);                      % Temperature
        V3 = R_I*T3/P3;                                                        % Spec. Volume
        a3 = sqrt(gamma_I*R_I*T3);                                             % Speed of sound
        
        delta_2 = delta_t - delta_1;                                           % Deflection angle
        omega_r = pi-asin(M2n_prime/M2);                                       % Angle of the reflected shock
        t_2 = y/u2 * sin(omega_r-omega_i+delta_1)/(sin(omega_r)*sin(omega_i)); % Time spent in zone 2
        
        t_3 = t_1;                                                             % Time spent in zone 3, arbitrary set to be equal to t_1
        
        if plt
            
            % Legend of the graph
            lgd = strcat('Mach of incident flow : ', num2str(M1), '; Particle x = ', num2str(x*100), 'cm ; y = ', num2str(y*100), 'cm');
            %ttl = lgd;
            %ttl = strcat('Mach of incident flow : ', num2str(M1), "; All particles");
            ttl = strcat('All Mach numbers; Particle x = ', num2str(x*100), 'cm ; y = ', num2str(y*100), 'cm');

            % Color of the graph
            col = varargin{2};
            
            t_plot = [0, t_1, t_1, t_1+t_2, t_1+t_2, t_1+t_2+t_3];             % Time array
            P_plot = [P1, P1, P2, P2, P3, P3];                                 % Pressure array
            T_plot = [T1, T1, T2, T2, T3, T3];                                 % Temperature array
            V_plot = [V1, V1, V2, V2, V3, V3];                                 % Spec volume array
            
            figure(1)
            subplot(3,1,1)
            pressure = plot(t_plot, P_plot);
            title({'Pressure evolution';ttl})
            xlabel("Time (s)")
            ylabel("Pressure (Pa)")
            pressure.Color = col;
            grid on
            grid minor
            hold on
            
            subplot(3,1,2)
            temperature = plot(t_plot, T_plot);
            title('Temperature evolution')
            xlabel("Time (s)")
            ylabel("Temperature (K)")
            temperature.Color = col;
            grid on
            grid minor
            hold on
            
            subplot(3,1,3)
            volume = plot(t_plot, V_plot, 'DisplayName', lgd);
            title('Specific volume evolution')
            xlabel("Time (s)")
            ylabel("Specific volume")
            volume.Color = col;
            legend();
            grid on
            grid minor
            hold on
            
            figure(2)
            [tmp]=RRRscheme(x,y,omega_i, omega_t, delta_1, delta_t, omega_r, col);
            title({'Trajectory of Lagrangian particle';ttl})
            hold on
        end
        solution = [t_1 t_2 t_3 0 0 ; P1 P2 P3 P4 P5 ; T1 T2 T3 T4 T5 ; V1 V2 V3 V4 V5];
    end
end
end