% This programm computes the evolution of pressure, temperature and
% specific volume of a Lagrangian particle when it travels through a RRE
% structure

% Inputs : mode = 1 : H2_O2 // He interface
%          x, y coordinates of the lag. particle at the beginning of the
%          simulation
%          mach : mach number of the incident shock
% Outputs : t_plot : time array
%          P, T, V : pressure temperature and specifoc volume

function [t, Ve, t_1, t_2, t_3, P1, T1]=RREevolution(mode, x, y, mach, plot)

if mode ~= 1
    disp("Unexisting mode")
else
    run H2_O2_conds.m
    r_0 = y - x/tan(omega_rad);
    if r_0 < 0
        disp("This code only follows particles in phase I!")
    else
        %% Constants before incident shock
        
        % Pressure P0; Temperature T0_I; Spec. Volume Vm0_I
        % --> set by H2_O2_conds.m
        
        M0 = mach; % Mach number
        M0_n = M0*sin(omega_rad); % Normal mach number
        
        t_0 = (experiment_dim+x)/(M0*a0_I); %time spent before shock
        % deflection set to be void
        
        %% Constants between incident shock and expansion fan
        pressJump = (1-gamma_I+2*gamma_I*M0_n^2)/(1+gamma_I);%pressure jump
        tempJump = (2+M0_n^2*(gamma_I-1))*(1-gamma_I+2*gamma_I*M0_n^2)/((gamma_I+1)^2*M0_n^2);%temperature jump
        volJump = (gamma_I-1)/(gamma_I+1) + 2/((gamma_I+1)*M0_n^2);%spec volume jump
        xi = pressJump;
        
        P1 = P0*pressJump;
        T1 = T0_I*tempJump;
        Vm1 = Vm0_I*volJump;
        
        M1_n = sqrt( (2 + (gamma_I - 1)*M0_n^2 / (1 - gamma_I + 2*gamma_I*M0_n^2) )); % Normal mach number
        % deflection angle
        d1 = atan( 2*cot(omega_rad)*(M0_n^2-1)/(M0^2*(gamma_I+cos(2*omega_rad))+2) );
        M1 = M1_n / sin(omega_rad-d1); % Mach number behind shock
        
        mu1=asin(1/M1); %angle of forward Mach wave
        
        %time spent between incident shock and expansion fan
        t_1 = r_0/(M1*a0_I)*( sin(omega_rad+d1)/tan(mu1) + cos(omega_rad+d1));
        
        %% Constants behind expansion fan
        
        % Prandtl-Meyer theory allows to calculate P2, T2, M2, V2 given the
        % total deflection behind the expansion (known thanks to transmitted
        % shock and shock polars)
        ny = 500; % number of points used to calculate polars
        %[i_xi, i_delta]=getPolar(M0, gamma_I, ny, 0);
        %semilogy(i_delta,i_xi,'b');
        %[out1,out2,out3]=getPolarPoint(i_xi, i_delta,0, xi);
        %graph_d1 = out2(2);
        %hold on;
        Mt = a0_I/a0_II*M0; % Transmitted mach number
        [t_xi, t_delta] = getPolar(Mt, gamma_II,ny,0); % Transmitted shock polar
        [r_xi, r_delta] = getExpPolar(M1, gamma_I, ny, xi, d1); % Expansion polar
        
        %semilogy(t_delta,t_xi,'r');
        %semilogy(r_delta,r_xi,'y');
        
        %[graph_xi, graph_delta] = getExpPolar(M1, gamma_I, ny, xi, graph_d1);
        % Intersection point between expansion polar and transmitted shock polar
        % allows to determine transmitted pressure jump and trabsmitted deviation
        [dev_t, xt] = getCrossPoint(t_delta, t_xi,r_delta, r_xi);
        
        % Deviation after expansion
        d2 = dev_t - d1;
        [~,nu_M1,~] = flowprandtlmeyer(gamma_I,M1);
        nu_M2 = d2 + nu_M1;
        % Mach number after expansion
        [M2, ~, ~] = flowprandtlmeyer(gamma_I, nu_M2, 'nu');
        
        %Pressure, temperature and spec. volume after expansion
        P2 = P0*xt;
        T2 = T1*(P2/P1)^((gamma_I-1)/gamma_I);
        V2 = R_I*T2/P2;
        
        %% Evolution in the expansion fan
        
        n_exp=1000; %number of points to solve differential equation
        ddelta=d2/(n_exp); %deviation step
        
        u1=M1*sqrt(gamma_I*R_I*T1); %speed behind expansion fan
        u1nexp=u1/M1; %speed normal to the beginning of the expansion fan
        u1texp=sqrt(u1^2-u1nexp^2); %speed tangent to the beginning of expansion fan
        t=zeros(1,n_exp); %initializing time array
        r=r_0*M1*sin(omega_rad-d1); %initializing radius i.e. distance to...
        %... point of incidence
        
        Msqexp=M1^2; %initializing square of Mach in expansion fan
        Te=zeros(1,n_exp); %initializing array of temperatures in expansion fan
        Pe=zeros(1,n_exp); %initializing array of pressures in expansion fan
        Ve=zeros(1,n_exp); %initializing array of spec. vols. in expansion fan
        Pexp=P1; %pressure at current position
        Texp=T1; %temperature at current position
        Vmexp=Vm1; %spec. vol. at current position
        
        for i=1:n_exp
            dMsqexp=2*Msqexp*(1+(gamma_I-1)/2*Msqexp)/sqrt(Msqexp-1)*ddelta; %square ...
            %... Mach variation
            
            %Computing pressure and temperature
            dPexp=-Pexp*(gamma_I/2)/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %pressure variation
            Pexp=Pexp+dPexp; %new pressure
            dTexp=-Texp*((gamma_I-1)/2)/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %temperature variation
            Texp=Texp+dTexp; %new temperature
            dVmexp=Vmexp/2*1/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %spec. vol. variaton
            Vmexp=Vmexp+dVmexp;
            Pe(i)=Pexp; %adding pressure to array
            Te(i)=Texp; %adding temperature to array
            Ve(i)=Vmexp; %adding spec. vol. to array
            
            %geometry calculations
            dtheta=ddelta+1/(2*sqrt(Msqexp^2-Msqexp))*dMsqexp; %angle variation...
            %... in cylindrical coordinates
            dt=r/u1nexp*dtheta; %time interval
            dr=u1texp*dt; %radius variation
            du1nexp=u1nexp/2*1/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %variation of...
            %... speed normal to small disturbance
            
            r=r+dr; %new radius
            t(i)=dt; %adding time to array
            u1nexp=u1nexp+du1nexp; %new normal speed
            
            %changing Msqexp
            Msqexp=Msqexp+dMsqexp;
        end
        
        %summing time intervals
        for i=2:n_exp
            t(i)=t(i-1)+t(i);
        end
        
        t_2 = t(end); %time spent in the expansion
        t_3 = 2e-6; %time spent after the expansion, arbitrary set to be 2ms
        
        if plot
            % Time array
            t_plot = zeros(1,1+2+n_exp+1);
            t_plot(1) = 0; % 1 initial point
            t_plot(2:3) = [t_0, t_0]; %2 jump points
            t_plot(4:3+n_exp) = (t_0+t_1)*ones(1,n_exp)+t; % n_exp expansion points
            t_plot(n_exp+4) = t_0+t_1+t_2+t_3; % 1 final point
            
            % Pressure array
            P_plot = zeros(1,1+2+n_exp+1);
            P_plot(1) = P0; % 1 initial point
            P_plot(2:3) = [P0, P1]; %2 jump points
            P_plot(4:3+n_exp) = Pe; % n_exp expansion points
            P_plot(n_exp+4) = P2; % 1 final point
            
            % Temperature array
            T_plot = zeros(1,1+2+n_exp+1);
            T_plot(1) = T0_I; % 1 initial point
            T_plot(2:3) = [T0_I, T1]; %2 jump points
            T_plot(4:3+n_exp) = Te; % n_exp expansion points
            T_plot(n_exp+4) = T2; % 1 final point
            
            % Spec volume array
            V_plot = zeros(1,1+2+n_exp+1);
            V_plot(1) = Vm0_I; % 1 initial point
            V_plot(2:3) = [Vm0_I, Vm1]; %2 jump points
            V_plot(4:3+n_exp) = Ve; % n_exp expansion points
            V_plot(n_exp+4) = V2; % 1 final point
            
            subplot(3,1,1)
            plot(t_plot, P_plot, '-.r')
            title('Pressure evolution')
            xlabel("Time (s)")
            ylabel("Pressure (Pa)")
            
            
            subplot(3,1,2)
            plot(t_plot, T_plot, '-.r')
            title('Temperature evolution')
            xlabel("Time (s)")
            ylabel("Temperature (K)")
            
            
            subplot(3,1,3)
            plot(t_plot, V_plot, '-.r')
            title('Specific volume evolution')
            xlabel("Time (s)")
            ylabel("Specific volume")
        end
    end
end