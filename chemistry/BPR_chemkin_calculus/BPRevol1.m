%Computation of the evolution of Lagragian particle of coordinates x and y
%in the coordinate system set in theory (see report)
%for a BPR refraction structure

%LP travels through incident shock, reflected shock and secondly refracted
%expansion
%at the beginning, pressure behind reflected shock is unknown : dichotomy
%is used to determine the unique deflection that allows existence of a
%shock, an expansion and P3=P4

function [solution, expansion] = BPRevol1 (cst, x, y, mach, omega_deg, press, temp_I, temp_II, varargin)

solution = -1;
expansion = -1;
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
        
        xi_t = machtoPressJump(M5n, gamma_II, pi/2);                         % Pressure jump
        P4 = xi_t*P5;                                                          % Pressure
        V4 = V5*machtoVolJump(M5n, gamma_II, pi/2);                           % Spec. volume
        T4 = T5*machtoTempJump(M5n, gamma_II, pi/2);                          % Temperature
        a4 = sqrt(gamma_II*R_II*T4);                                             % Speed of sound
        delta_t = postShockDeflection(M5, gamma_II, omega_t);                   % Deflection angle behind the shock
        
        %M4 = M4n / sin(omega_t - delta_t);                                     % Mach number
        %u4 = M4*a4;
        
        %% Zone 2b and 3 to be determined by testing different angles omega_r
        % delta_2 + delta_3 = delta_t - delta_1
        
        % Test for omega_r in [0 pi/2]
        n_om = 10000;
        omega_r_array = linspace(0.1, pi/10,n_om);
        
        delta_2_array = zeros(1,n_om);
        delta_3_array = zeros(1,n_om);
        M2an_prime_array = zeros(1,n_om);
        M2bn_array = zeros(1,n_om);
        M2b_array = zeros(1,n_om);
        P2b_array = zeros(1,n_om);
        M3_array = zeros(1,n_om);
        P3_array = zeros(1,n_om);
        
        for i = 1:1:n_om
            delta_2_array(i) = postShockDeflection(M2a, gamma_I, omega_r_array(i));
            delta_3_array(i) = delta_t - delta_1 - delta_2_array(i);
            
            M2an_prime_array(i) = M2a*sin(omega_r_array(i));
            M2bn_array(i) = machtoPSNormalMach(M2an_prime_array(i), gamma_I, pi/2);
            M2b_array(i) = M2bn_array(i)/sin(omega_r_array(i) - delta_2_array(i));
            P2b_array(i) = P2a*machtoPressJump(M2an_prime_array(i), gamma_I, pi/2);
        
            [~,nu_M2b, ~] = flowprandtlmeyer(gamma_I, M2b_array(i), 'mach');
            nu_M3 = delta_3_array(i) + nu_M2b;
            [M3_array(i), ~, ~] = flowprandtlmeyer(gamma_I, nu_M3, 'nu');
        
            P3_array(i) = P2b_array(i)*pressureDropExp(M2b_array(i), M3_array(i), gamma_I);
        end
        [diff, index] = min(abs(P3_array - P4));
        
        if diff/P4 < 1e-3
            omega_r = omega_r_array(index);
            delta_2 = delta_2_array(index);
            delta_3 = delta_3_array(index);
            M2b = M2b_array(index);
            P2b = P2b_array(index);
            M3 = M3_array(index);
            P3 = P3_array(index);
        else
            disp('Solution not found for omega_r < pi/2')
            %return;
        end
        
        %% Zone 2b
        % delta_2 and omega_r = already known
        % M2bn and M2b = already known
        
        t2a = y/u2a * sin(omega_r-omega_i+delta_1)/(sin(omega_r)*sin(omega_i)); % Time spent in zone 2a
        
        xi_r = P2b/P2a;
        V2b = V2a*machtoVolJump(M2a, gamma_I, pi-omega_r);
        T2b = T2a*machtoTempJump(M2a, gamma_I, pi-omega_r);
        a2b = sqrt(gamma_I*R_I*T2b);
        u2b = M2b*a2b;
        mu_2b = asin(1/M2b);
        r=y*sin(omega_i-delta_1)/(sin(omega_r)*sin(omega_i));                    % Position of particle on reflected shock
        
        t2b = r/u2b * sin(mu_2b + omega_r - delta_2)/sin(mu_2b);
        r2 = r * sin(omega_r - delta_2)/sin(mu_2b);
        %% Zone 3
        % delta_3 and M3 = already known
        
        T3 = T2b*temperatureDropExp(M2b, M3, gamma_I);
        V3 = V2b*T3*P2b/(T2b*P3);
        a3 = sqrt(gamma_I*R_I*T3);
        
        t3 = t1; % t3 arbitrary set to be equal to t1
        %% Zone e (isentropic evolution)
        
        n_exp=1000;                                                            % Number of points to solve differential equation
        ddelta=delta_3/(n_exp);                                                % deviation step
        u2bnexp=u2b/M2b;                                                          % Normal velocity at the beginning of the expansion fan
        u2btexp=sqrt(u2b^2-u2bnexp^2);                                            % Tangential velocity at the beginning of expansion fan
        te=zeros(1,n_exp);                                                     % Time array init.
        
        Me2=M2b^2;                                                              % Square of Mach number init.
        Te=zeros(1,n_exp);                                                     % Temperature array init.
        Pe=zeros(1,n_exp);                                                     % Pressure array init.
        Ve=zeros(1,n_exp);                                                     % Spec. vol. array init.
        Pexp=P2b;                                                               % Pressure at current position
        Texp=T2b;                                                               % Temperature at current position
        Vexp=V2b;                                                               % Spec. vol. at current position
        
        Pe(1)=Pexp;
        Te(1)=Texp;
        Ve(1)=Vexp;
        
        for i=2:n_exp
            dMsqexp=2*Me2*(1+(gamma_I-1)/2*Me2)/sqrt(Me2-1)*ddelta;            % Square Mach variation
            
            % Computing pressure and temperature
            dPexp=pressVarExp(Pexp, Me2, dMsqexp, gamma_I);                    % Pressure variation
            Pexp=Pexp+dPexp;                                                   % New pressure
            dTexp=tempVarExp(Texp, Me2, dMsqexp, gamma_I);                     % Temperature variation
            Texp=Texp+dTexp;                                                   % New temperature
            dVmexp=volVarExp(Vexp, Me2, dMsqexp, gamma_I);                     % Spec. vol. variaton
            Vexp=Vexp+dVmexp;                                                  % New specific volume
            Pe(i)=Pexp;                                                        % Adding pressure to array
            Te(i)=Texp;                                                        % Adding temperature to array
            Ve(i)=Vexp;                                                        % Adding spec. vol. to array
            
            % Geometry calculations
            dtheta=ddelta+1/(2*sqrt(Me2^2-Me2))*dMsqexp;                       % Angle variation in cylindrical coordinates
            dt=r2/u2bnexp*dtheta;                                                % time interval
            dr2=u2btexp*dt;                                                      % Radius variation
            du1nexp=u2bnexp/2*1/(1+(gamma_I-1)/2*Me2)*dMsqexp;                  % Variation of velocity normal to small disturbance
            
            r2=r2+dr2;                                                            % New radius
            te(i)=te(i-1)+dt;                                                  % Adding time to array
            u2bnexp=u2bnexp+du1nexp;                                             % New normal speed
            
            Me2=Me2+dMsqexp;                                                   % Changing Msqexp
        end
        
        mu_3 = asin(1/sqrt(Me2));
        texp = te(end); % Time spent in the expansion
        
        if plt
            % Legend of the graph
            lgd = strcat('Mach of incident flow : ', num2str(M1), '; Particle x = ', num2str(x*100), 'cm ; y = ', num2str(y*100), 'cm');
            %ttl = lgd;
            ttl = strcat('Mach of incident flow : ', num2str(M1), "; All particles");
            %ttl = strcat('All Mach numbers; Particle x = ', num2str(x*100), 'cm ; y = ', num2str(y*100), 'cm');
            
            % Color of the graph
            col = varargin{2};
            
            % Time array
            t_plot = zeros(1,1+2+2+n_exp+1);
            t_plot(1) = 0; % 1 initial point
            t_plot(2:3) = [t1, t1]; %2 jump points
            t_plot(4:5) = [t1+t2a, t1+t2a]; %2 jump points
            t_plot(5:4+n_exp) = (t1+t2a+t2b)*ones(1,n_exp)+te; % n_exp expansion points
            t_plot(n_exp+5) = t1+t2a+t2b+texp+t3; % 1 final point
            
            % Pressure array
            P_plot = zeros(1,1+2+2+n_exp+1);
            P_plot(1) = P1; % 1 initial point
            P_plot(2:3) = [P1, P2a]; %2 jump points
            P_plot(4:5) = [P2a, P2b]; %2 jump points
            P_plot(5:4+n_exp) = Pe; % n_exp expansion points
            P_plot(n_exp+5) = P3; % 1 final point
            
            % Temperature array
            T_plot = zeros(1,1+2+2+n_exp+1);
            T_plot(1) = T1; % 1 initial point
            T_plot(2:3) = [T1, T2a]; %2 jump points
            T_plot(4:5) = [T2a, T2b]; %2 jump points
            T_plot(5:4+n_exp) = Te; % n_exp expansion points
            T_plot(n_exp+5) = T3; % 1 final point
            
            % Spec volume array
            V_plot = zeros(1,1+2+2+n_exp+1);
            V_plot(1) = V1; % 1 initial point
            V_plot(2:3) = [V1, V2a]; %2 jump points
            V_plot(4:5) = [V2a, V2b]; %2 jump points
            V_plot(5:4+n_exp) = Ve; % n_exp expansion points
            V_plot(n_exp+5) = V3; % 1 final point
            
            figure(1)
            subplot(3,1,1)
            pressure = plot(t_plot, P_plot);
            title({ttl; 'Pressure evolution'})
            xlabel("Time (s)")
            ylabel("Pressure (Pa)")
            pressure.Color = col;
            hold on
            grid on
            grid minor
            
            
            
            subplot(3,1,2)
            temperature = plot(t_plot, T_plot);
            title('Temperature evolution')
            xlabel("Time (s)")
            ylabel("Temperature (K)")
            temperature.Color = col;
            hold on
            grid on
            grid minor
            
            
            subplot(3,1,3)
            volume = plot(t_plot, V_plot, 'DisplayName', lgd);
            title('Specific volume evolution')
            xlabel("Time (s)")
            ylabel("Specific volume")
            volume.Color = col;
            legend();
            hold on
            grid on
            grid minor
        end
        
        solution = [t1 t2a t2b t3 0 0 ; P1 P2a P2b P3 P4 P5 ; T1 T2a T2b T3 T4 T5 ; V1 V2a V2b V3 V4 V5];
        expansion = [te ; Pe ; Te ; Ve];
    end
end