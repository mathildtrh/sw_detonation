% This programm computes the evolution of the specific volume of a
% lagrangian particle x,y when it travels through a RRE struture

% mode allows you to select which reactive-inert interface you want to
% simulate
% mach allows you to choose the mach number of the incident flow

% t = time array in the expansion
% V = specific volume array in the expansion
% P1, T1 pressure and temperature behind the incident shock

function [ts, Vms, t_1, t_2, t_3, P1, T1]=partThroughRRE(mode, x, y, mach)

% Selection of the reactive-inert interface
if mode == 1
    H2_O2_conds;
else
    "Mode inexistant"
end

%distance to point of incidence, when shock reaches particle
r_0=y-x/tan(omega_rad); %has to be positive!
if r_0<0
    "This code only follows particles in phase I!"
else
    %% Calculation of flow behind incident shock
    
    % Normal square mach before shock
    Msq0_I_N = mach^2*sin(pi - omega_rad)^2;
    
    % Pressure jump
    xi = machToXi(mach,gamma_I, pi - omega_rad);
    P1 = P0*xi/101325; %conversion in atm
    
    % Deflection angle behind the incident shock
    delta_1=atan(sqrt(tanDefSq(xi,Msq0_I_N,gamma_I)));
    
    % Temperature jump
    T1 = T0_I*xiToTempJ(xi, gamma_I); %already in K
    
    % Specific volume jump
    Vm1 = Vm0_I*xiToVolJump(xi, gamma_I);
    
    % Normal square mach number behind the incident shock
    Msq1N = ( 2 + (gamma_I-1)*Msq0_I_N ) / (1 - gamma_I + 2*gamma_I*Msq0_I_N);
    
    % Square mach number behind the incident shock
    Msq1 = Msq1N/sin(pi - omega_rad - delta_1);
    
    delta_1=atan(sqrt(tanDefSq(xi,Msq0_I_N,gamma_I)));
    mu_1=asin(1/Msq0_I_N); %velocity vector-expansion fan angle at beginning...
    %... of expansion fan
    t_1=1/sqrt(gamma_I*R_I*T0_I)*r_0*sin(omega_rad-delta_1+mu_1); %time spent...
    %between the incident shock and the expansion
    
    
    %% When in zone 2, behind the expansion fan
    Mr=sqrt(Msq1);%Mach when arriving at reflection
    
    %Computing expansion deflection and pressure jump
    Mt=sqrt((gamma_I*T0_I*mu_II)/(gamma_II*T0_II*mu_I))*Msq0_I_N;
    npts=200; %number of points
    xis=logspace(0,log(xi)/log(10),npts); %pressure jumps when going through...
    %... expansion fan
    deltas_r=ones(1,npts)*delta_1; %initializing array of expansion deflections
    deltas_t=zeros(1,npts); %initializing array of transmitted wave deflections
    for i=1:npts
        M_post_exp=sqrt(2/(gamma_I-1)*((1+(gamma_I-1)/2*Mr^2)*(xis(i)/xi)^...
            ((1-gamma_I)/gamma_I)-1)); %post-expansion Mach at given pressure jump
        deltas_r(i)= deltas_r(i) + prandtlMeyer(M_post_exp,gamma_I)-...
            prandtlMeyer(Mr,gamma_I); %computing deltas_r
        
        deltas_t(i)=atan(sqrt(tanDefSq(xis(i),Mt,gamma_II))); %computing...
        %... transmitted wave deflection
    end
    
    def_diff=abs(deltas_r-deltas_t);
    [m,k]=min(def_diff);
    delta_2=deltas_r(k)-delta_1; %deviation at end of expansion
    xi_t=xis(k); %pressure jump after expansion
    M_post_exp=sqrt(2/(gamma_I-1)*((1+(gamma_I-1)/2*Mr^2)*(xi_t/xi)^...
        ((1-gamma_I)/gamma_I)-1)); %Mach after expansion
    exp_temp_jump=(1+(gamma_I-1)/2*Mr^2)...
        /(1+(gamma_I-1)/2*M_post_exp^2); %temperature jump after expansion
    T_post_exp=T1*exp_temp_jump; %temperature after expansion
    exp_spec_vol_jump=((1+(gamma_I-1)/2*Mr^2)...
        /(1+(gamma_I-1)/2*M_post_exp^2))^(1/(1-gamma_I)); %spec. vol. jump...
    %... after expansion
    Vm_post_exp=Vm1*exp_spec_vol_jump;
    
    %% In expansion fan
    no_dev_pts=1000; %number of points to solve differential equation
    ddelta=delta_2/(no_dev_pts); %deviation step

    u1=Mr*sqrt(gamma_I*R_I*T1); %speed before expansion fan
    u1nexp=u1/Mr; %speed normal to the beginning of the expansion fan
    u1texp=sqrt(u1^2-u1nexp^2); %speed tangent to the beginning of expansion fan
    tsexp=zeros(1,no_dev_pts); %initializing time array
    r=r_0*Mr*sin(omega_rad-delta_1); %initializing radius i.e. distance to...
    %... point of incidence
    
    Msqexp=Mr^2; %initializing square of Mach in expansion fan
    Texps=zeros(1,no_dev_pts); %initializing array of temperatures in expansion fan
    Pexps=zeros(1,no_dev_pts); %initializing array of pressures in expansion fan
    Vmexps=zeros(1,no_dev_pts); %initializing array of spec. vols. in expansion fan
    Pexp=P0*xi; %pressure at current position
    Texp=T1; %temperature at current position
    
    for i=1:no_dev_pts
        dMsqexp=2*Msqexp*(1+(gamma_I-1)/2*Msqexp)/sqrt(Msqexp-1)*ddelta; %square ...
        %... Mach variation
        
        %Computing pressure and temperature
        dPexp=-Pexp*(gamma_I/2)/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %pressure variation
        Pexp=Pexp+dPexp; %new pressure
        dTexp=-Texp*((gamma_I-1)/2)/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %temperature variation
        Texp=Texp+dTexp; %new temperature
        Vmexp=R_I*Texp/Pexp; %new spec. volume, with perfect gases law
        
        Pexps(i)=Pexp; %adding pressure to array
        Texps(i)=Texp; %adding temperature to array
        Vmexps(i)=Vmexp; %adding spec. vol. to array
        
        %geometry calculations
        dtheta=ddelta+1/(2*sqrt(Msqexp^2-Msqexp))*dMsqexp; %angle variation...
        %... in cylindrical coordinates
        dt=r/u1nexp*dtheta; %time interval
        dr=u1texp*dt; %radius variation
        du1nexp=u1nexp/2*1/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %variation of...
        %... speed normal to small disturbance
        
        r=r+dr; %new radius
        tsexp(i)=dt; %adding time to array
        u1nexp=u1nexp+du1nexp; %new normal speed
        
        %changing Msqexp
        Msqexp=Msqexp+dMsqexp;
    end
    
    %summing time intervals
    for i=2:no_dev_pts
        tsexp(i)=tsexp(i-1)+tsexp(i);
    end
    
    t_2 = tsexp(end);
    t_3 = t_1;
    
    %construction of time array
    %Chemkin needs time to start right before the expansion  
    ts=zeros(1,no_dev_pts+1);
    ts(1) = 0;
    ts(2:no_dev_pts+1)=tsexp;
    
    % construction of specific volume array
    Vms=zeros(1,no_dev_pts+1);
    Vms(1)=Vm1;
    Vms(2:no_dev_pts+1)= Vmexps;
    
    
    if ~(length(ts) == length(Vms))
        "PB DE DIMENSION"
    end
end

end