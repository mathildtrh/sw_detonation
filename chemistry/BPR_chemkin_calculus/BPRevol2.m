%Computation of the evolution of Lagragian particle of coordinates x and y
%in the coordinate system set in theory (see report)
%for a BPR refraction structure

%LP travels through incident shock, reflected shock and secondly refracted
%expansion --> this time, let's try without refracted expansion

%at the beginning, pressure behind reflected shock is unknown : dichotomy
%is used to determine the unique deflection that allows existence of a
%shock and P2b=P4

function [solution] = BPRevol2 (cst, x, y, mach, omega_deg, press, temp_I, temp_II, varargin)

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
        
        ui = mach*a1;
        M1 = mach/sin(omega_i);                                                % Mach number in zone 1
        u1 = M1*a1;
        t1 = 1/u1 * (x - y/tan(omega_i));                                      % Time spent in zone 1
        
        %% Zone 5
        
        P5 = press;                                                            % Pressure
        T5 = temp_II;                                                          % Temperature
        V5 = R_II*T5/P5;                                                       % Spec. volume
        a5 = sqrt(gamma_II*R_II*T5);                                           % Speed of sound
        
        %% Zone 2a (oblique shock from 1)
        
        M2an = machtoPSNormalMach(M1, gamma_I, omega_i);                        % Mach number (2) normal to incident shock
        xi_i = machtoPressJump(M1, gamma_I, omega_i);                          % Pressure jump
        P2a = xi_i*P1;                                                          % Pressure
        V2a = V1*machtoVolJump(M1, gamma_I, omega_i);                           % Spec. volume
        T2a = T1*machtoTempJump(M1, gamma_I, omega_i);                          % Temperature
        a2a = sqrt(gamma_I*R_I*T2a);                                             % Speed of sound
        delta_1 = postShockDeflection(M1, gamma_I, omega_i);                   % Deflection angle behind the shock
        
        M2a = M2an / sin(omega_i - delta_1);                                     % Mach number
        u2a = M2a*a2a;
        
        %% Zone 4 (oblique shock from 5 + piston theory)
        b = (gamma_II+1)/(gamma_I+1)*(ui^2-a1^2)/ui;                           % Temporary value
        ut = 0.5*(b+sqrt(b^2+4*a5^2));                                         % Normal speed of transmitted shock
        omega_t = asin(sin(omega_i)*ut/ui);                                    % Angle of transmitted shock
        M5n = ut/a5;                                                           % Normal Mach number
        M5 = M5n/sin(omega_t);                                                 % Mach number
        
        xi_t = machtoPressJump(M5, gamma_II, omega_t);                         % Pressure jump
        P4 = xi_t*P5;                                                          % Pressure
        V4 = V5*machtoVolJump(M5, gamma_II, omega_t);                           % Spec. volume
        T4 = T5*machtoTempJump(M5, gamma_II, omega_t);                          % Temperature
        a4 = sqrt(gamma_II*R_II*T4);                                             % Speed of sound
        delta_t = postShockDeflection(M5, gamma_II, omega_t);                   % Deflection angle behind the shock
        
        %M4 = M4n / sin(omega_t - delta_t);                                     % Mach number
        %u4 = M4*a4;
        
        %% Let's suppose that zone 3 doesn't exist
        %% Zone 2b to be determined by solving the angle of shock
        % delta_2 = delta_t - delta_1
        
        P2b = P4;
        xi_r = P2b/P2a;
        M2an_prime = sqrt( xiToSqMach(xi_r, gamma_I, pi/2) );
        M2bn = machtoPSNormalMach(M2an_prime, gamma_I, pi/2);
        V2b = V2a*machtoVolJump(M2an_prime, gamma_I, pi/2);
        T2b = T2a*machtoTempJump(M2an_prime, gamma_I, pi/2);
        a2b = sqrt(gamma_I*R_I*T2b);
        
        omega_r = pi-asin(M2an_prime/M2a);
        delta_2 = postShockDeflection(M2a, gamma_I, pi-omega_r);
        
        if abs(delta_2 + delta_1-delta_t)/delta_t > 1e-3
            disp('Discrepancies on deflection angles')
        end
        
        M2b = M2bn/sin(omega_r - delta_2);
        t2a = y/u2a * sin(omega_r-omega_i+delta_1)/(sin(omega_r)*sin(omega_i)); % Time spent in zone 2a

        r=y*sin(omega_i-delta_1)/(sin(omega_r)*sin(omega_i));                   % Position of particle on reflected shock
        t2b = t1; %arbitrary set

        if plt
            
            % Legend of the graph
            lgd = strcat('Mach of incident flow : ', num2str(M1), '; Particle x = ', num2str(x*100), 'cm ; y = ', num2str(y*100), 'cm');
            %ttl = lgd;
            %ttl = strcat('Mach of incident flow : ', num2str(M1), "; All particles");
            ttl = strcat('All Mach numbers; Particle x = ', num2str(x*100), 'cm ; y = ', num2str(y*100), 'cm');

            % Color of the graph
            col = varargin{2};
            
            t_plot = [0, t_1, t_1, t_1+t_2a, t_1+t_2a, t_1+t_2a+t_2b];             % Time array
            P_plot = [P1, P1, P2a, P2a, P2b, P2b];                                 % Pressure array
            T_plot = [T1, T1, T2a, T2a, T2b, T2b];                                 % Temperature array
            V_plot = [V1, V1, V2a, V2a, V2b, V2b];                                 % Spec volume array
            
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
        end
        
        solution = [t1 t2a t2b 0 0 ; P1 P2a P2b P4 P5 ; T1 T2a T2b T4 T5 ; V1 V2a V2b V4 V5];
    end
end