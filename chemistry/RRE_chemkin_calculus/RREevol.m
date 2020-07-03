%Computation of the evolution of Lagragian particle of coordinates x and y
%in the coordinate system set in theory (see report)

function [solution, expansion] = RREevol (cst, x, y, mach, omega_deg, press, temp_I, temp_II, plt)

solution = -1;
expansion = -1;
gamma_I = cst(1);
gamma_II = cst(2);
R_I = cst(3);
R_II = cst(4);

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
        
        ui = mach*a1;                                                          % Speed of the shock
        u1 = sqrt(ui^2*(1+(tan(pi-omega_i))^2));                               % u1 = -u_node -> change of frame of reference
        M1 = u1/a1;                                                            % Mach number in zone 1
        
        t1 = 1/u1 * (x - y/tan(omega_i));                                      % Time spent in zone 1
        
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
        mu_2 = asin(1/M2);                                                     % Forward Mach wave angle
        t2 = y/u2 * sin(mu_2+omega_i-delta_1)/(sin(mu_2)*sin(omega_i));        % Time spent in zone 2
        
        %% Zone 4 (oblique shock from 5)
        
        ny= 500;                                                               % Number of points of computation for polar
        [refl_x, refl_d] = getExpPolar(M2, gamma_I, ny, xi_i, delta_1);        % Polar of reflected expansion
        [trans_x, trans_d] = getPolar(M5, gamma_II, ny);                       % Polar of transmitted shock
        
        [delta_t, xi_t] = getCrossPoint(trans_d, trans_x, refl_d, refl_x, 1);  % Get the intersection point between the two - see Yann's report
        
        P4 = xi_t*P5;                                                          % Pressure
        V4 = V5*xiToVolJump(xi_t, gamma_II);                                   % Spec. Volume
        T4 = T5*xi_t*xiToVolJump(xi_t, gamma_II);                              % Temperature
        a4 = sqrt(gamma_II*R_II*T4);                                           % Speed of sound
        omega_t = asin(xiToSqMach(xi_t, gamma_II, pi/2)/M5); % omega_t = asin(M5n/M5)
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
        
        for i=1:n_exp
            dMsqexp=2*Me2*(1+(gamma_I-1)/2*Me2)/sqrt(Me2-1)*ddelta;            % Square Mach variation
            
            %Computing pressure and temperature
            dPexp=-Pexp*(gamma_I/2)/(1+(gamma_I-1)/2*Me2)*dMsqexp; %pressure variation
            Pexp=Pexp+dPexp; %new pressure
            dTexp=-Texp*((gamma_I-1)/2)/(1+(gamma_I-1)/2*Me2)*dMsqexp; %temperature variation
            Texp=Texp+dTexp; %new temperature
            dVmexp=Vexp/2*1/(1+(gamma_I-1)/2*Me2)*dMsqexp; %spec. vol. variaton
            Vexp=Vexp+dVmexp;
            Pe(i)=Pexp; %adding pressure to array
            Te(i)=Texp; %adding temperature to array
            Ve(i)=Vexp; %adding spec. vol. to array
            
            %geometry calculations
            dtheta=ddelta+1/(2*sqrt(Me2^2-Me2))*dMsqexp; %angle variation in cylindrical coordinates
            dt=r/u2nexp*dtheta; %time interval
            dr=u2texp*dt; %radius variation
            du1nexp=u2nexp/2*1/(1+(gamma_I-1)/2*Me2)*dMsqexp; %variation of velocity normal to small disturbance
            
            r=r+dr; %new radius
            te(i)=dt; %adding time to array
            u2nexp=u2nexp+du1nexp; %new normal speed
            
            %changing Msqexp
            Me2=Me2+dMsqexp;
        end
        
        %summing time intervals
        for i=2:n_exp
            te(i)=te(i-1)+te(i);
        end
        
        mu_3 = asin(1/sqrt(Me2));
        texp = te(end); %time spent in the expansion
        
        if plt
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
            plot(t_plot, P_plot, 'r')
            title('Pressure evolution')
            xlabel("Time (s)")
            ylabel("Pressure (Pa)")
            
            
            subplot(3,1,2)
            plot(t_plot, T_plot, 'r')
            title('Temperature evolution')
            xlabel("Time (s)")
            ylabel("Temperature (K)")
            
            
            subplot(3,1,3)
            plot(t_plot, V_plot, 'r')
            title('Specific volume evolution')
            xlabel("Time (s)")
            ylabel("Specific volume")
            
            figure(2)
            tmp = RREscheme(x,y,omega_i, omega_t, delta_1, delta_t, mu_2, mu_3);
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