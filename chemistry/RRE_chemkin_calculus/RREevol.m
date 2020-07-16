%Computation of the evolution of Lagragian particle of coordinates x and y
%in the coordinate system set in theory (see report)

function [solution, expansion] = RREevol (cst, x, y, mach, omega_deg, press, temp_I, temp_II, varargin)

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
        
        M1 = mach/sin(omega_i);                                                % Mach number in zone 1
        u1 = M1*a1;                                                            % u1 = -u_node -> change of frame of reference
        
        t1 = 1/u1 * (x - y/tan(omega_i));                                      % Time spent in zone 1
        
        %% Zone 5
        
        P5 = press;                                                            % Pressure
        T5 = temp_II;                                                          % Temperature
        V5 = R_II*T5/P5;                                                       % Spec. volume
        a5 = sqrt(gamma_II*R_II*T5);                                           % Speed of sound
        M5 = a1/a5*M1;                                                         % Mach number
        
        %% Zone 2 (oblique shock from 1)
        
        M2n = machtoPSNormalMach(M1, gamma_I, omega_i);                        % Mach number (2) normal to incident shock
        inc_x = machtoPressJump(M1, gamma_I, omega_i);                          % Pressure jump
        P2 = inc_x*P1;                                                          % Pressure
        V2 = V1*machtoVolJump(M1, gamma_I, omega_i);                           % Spec. volume
        T2 = T1*machtoTempJump(M1, gamma_I, omega_i);                          % Temperature
        a2 = sqrt(gamma_I*R_I*T2);                                             % Speed of sound
        delta_1 = postShockDeflection(M1, gamma_I, omega_i);                   % Deflection angle behind the shock
        
        M2 = M2n / sin(omega_i - delta_1);                                     % Mach number
        u2 = M2*a2;                                                            % Velocity of flow
        mu_2 = asin(1/M2);                                                     % Forward Mach wave angle
        t2 = y/u2 * sin(mu_2+omega_i-delta_1)/(sin(mu_2)*sin(omega_i));        % Time spent in zone 2
        
        %% Zone 4 (oblique shock from 5)
        
        ny= 500;                                                               % Number of points of computation for polar
        [refl_x, refl_d] = getExpPolar(M2, gamma_I, ny, inc_x, delta_1);       % Polar of reflected expansion
        [trans_x, trans_d] = getPolarMathilde(M5, gamma_II, ny);                       % Polar of transmitted shock

        [delta_t, xi_t] = getCrossPoint(trans_d, trans_x, refl_d, refl_x, 1);  % Get the intersection point between the two - see Yann's report
        
        P4 = xi_t*P5;                                                          % Pressure
        V4 = V5*xiToVolJump(xi_t, gamma_II);                                   % Spec. Volume
        T4 = T5*xi_t*xiToVolJump(xi_t, gamma_II);                              % Temperature
        a4 = sqrt(gamma_II*R_II*T4);                                           % Speed of sound
        omega_t = asin(sqrt(xiToSqMach(xi_t, gamma_II, pi/2))/M5);             % omega_t = asin(M5n/M5)
        
        % Plot polars if needed
        figure()
        [inc_x, inc_d]=getPolarMathilde(M1, gamma_I, ny, 0);
        inc = semilogy(inc_d*180/pi,inc_x,'r', 'DisplayName', 'Incident shock');
        hold on;
        trans = semilogy(trans_d*180/pi,trans_x,'g', 'DisplayName', 'Transmitted shock');
        refl = semilogy(refl_d*180/pi,refl_x,'y', 'DisplayName', 'Reflected expansion');
        inter = semilogy([delta_t*180/pi], [xi_t], '+k', 'DisplayName', 'Intersection point : (\delta_t, \xi_t)');
        inc.LineWidth = 2;
        trans.LineWidth = 2;
        refl.LineWidth = 2;
        inter.LineWidth = 2;
        grid on
        grid minor
        xlabel('\omega (in deg)')
        ylabel('\xi')
        legend('Location', 'eastoutside')
        title({'Polar figure for RRE'; strcat('M1 = ', num2str(M1), '; Particle x = ', num2str(x*100), 'cm')})
        
        %% Zone 3 (Prandtl-Meyer expansion from 2)
        
        delta_2 = delta_t - delta_1;                                           % Post expansion deflection
        [~, nu_M2, ~] = flowprandtlmeyer(gamma_I,M2);                          % Value of Prandtl-Meyer function for M2
        nu_M3 = delta_2 + nu_M2;                                               % Value of Prandtl-Meyer function for M3
        [M3, ~, ~] = flowprandtlmeyer(gamma_I, nu_M3, 'nu');                   % Mach number (3)
        
        P3 = P4;                                                               % Pressure
        T3 = T2*(P3/P2)^((gamma_I-1)/gamma_I);                                 % Temperature
        V3 = R_I*T3/P3;                                                        % Spec. Volume
        a3 = sqrt(gamma_I*R_I*T3);                                             % Speed of sound
        u3 = M3*a3;                                                            % Velocity of flow
        t3 = t1;                                                             % Time spent in zone 3 (arbitrary)
        
        %% Zone e (isentropic evolution)
        
        n_exp=1000;                                                            % Number of points to solve differential equation
        ddelta=delta_2/(n_exp);                                                % deviation step
        u2nexp=u2/M2;                                                          % Normal velocity at the beginning of the expansion fan
        u2texp=sqrt(u2^2-u2nexp^2);                                            % Tangential velocity at the beginning of expansion fan
        te=zeros(1,n_exp);                                                     % Time array init.
        r=y*sin(omega_i-delta_1)/(sin(mu_2)*sin(omega_i));                     % Radius init.
        
        Me2=M2^2;                                                              % Square of Mach number init.
        Te=zeros(1,n_exp);                                                     % Temperature array init.
        Pe=zeros(1,n_exp);                                                     % Pressure array init.
        Ve=zeros(1,n_exp);                                                     % Spec. vol. array init.
        Pexp=P2;                                                               % Pressure at current position
        Texp=T2;                                                               % Temperature at current position
        Vexp=V2;                                                               % Spec. vol. at current position
        
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
            dt=r/u2nexp*dtheta;                                                % time interval
            dr=u2texp*dt;                                                      % Radius variation
            du1nexp=u2nexp/2*1/(1+(gamma_I-1)/2*Me2)*dMsqexp;                  % Variation of velocity normal to small disturbance
            
            r=r+dr;                                                            % New radius
            te(i)=te(i-1)+dt;                                                  % Adding time to array
            u2nexp=u2nexp+du1nexp;                                             % New normal speed
            
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
            t_plot = zeros(1,1+2+n_exp+1);
            t_plot(1) = 0; % 1 initial point
            t_plot(2:3) = [t1, t1]; %2 jump points
            t_plot(4:3+n_exp) = (t1+t2)*ones(1,n_exp)+te; % n_exp expansion points
            t_plot(n_exp+4) = t1+t2+texp+t3; % 1 final point
            
            % Pressure array
            P_plot = zeros(1,1+2+n_exp+1);
            P_plot(1) = P1; % 1 initial point
            P_plot(2:3) = [P1, P2]; %2 jump points
            P_plot(4:3+n_exp) = Pe; % n_exp expansion points
            P_plot(n_exp+4) = P3; % 1 final point
            
            % Temperature array
            T_plot = zeros(1,1+2+n_exp+1);
            T_plot(1) = T1; % 1 initial point
            T_plot(2:3) = [T1, T2]; %2 jump points
            T_plot(4:3+n_exp) = Te; % n_exp expansion points
            T_plot(n_exp+4) = T3; % 1 final point
            
            % Spec volume array
            V_plot = zeros(1,1+2+n_exp+1);
            V_plot(1) = V1; % 1 initial point
            V_plot(2:3) = [V1, V2]; %2 jump points
            V_plot(4:3+n_exp) = Ve; % n_exp expansion points
            V_plot(n_exp+4) = V3; % 1 final point
            
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
            
            
            
            
            figure(2)
            tmp = RREscheme(x,y,omega_i, omega_t, delta_1, delta_t, mu_2, mu_3, col);
            title({'Trajectory of Lagrangian particle';ttl})
            hold on
            
        end

        %% Verification
        % false if the end of the calculated expansion is far from values computed in zone 3
        
        check = [false, false, false];
        if abs(Pe(end)-P3)/P3 < 1e-4
            check(1) = true;
        end
        if abs(Te(end)-T3)/T3 < 1e-4
            check(2) = true;
        end
        if abs(Ve(end)-V3)/V3 < 1e-4
            check(3) = true;
        end
        % disp(check);
        
        %% General Solution
        % each column of "solution" array is a zone
        % lines : 1-> Time / 2-> Pressure / 3-> Temperature / 4-> Spec. Volume
        
        % idem for "expansion" array, each column is an interval between two
        % Mach waves
        
        solution = [t1 t2 t3 0 0 ; P1 P2 P3 P4 P5 ; T1 T2 T3 T4 T5 ; V1 V2 V3 V4 V5];
        expansion = [te ; Pe ; Te ; Ve];
        
        % disp("xi_i ="); disp(xi_i); disp("xi_t="); disp(xi_t);
        % disp("delta_1 ="); disp(delta_1); disp("delta_t="); disp(delta_t);
        
    end
end
end