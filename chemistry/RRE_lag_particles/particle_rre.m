% Program inspired from particle_evol.m programm
% which can be found in examples/lagrangian_particle_rre
% to compute normalized specific volume evolution 
% of a lagrangian particle, depending on its initial position
% for a certain reactive-inert gas interface


% Inputs
    % mode = 1 --> shock from H2_O2 phase to He phase
    % mode = 2 --> shock from CO2 phase to He phase
    % x = x coordinate of the lagrangian particle at study
    % y = y coordinate of the lagrangian particle at study
% Outputs
    % ts = array of time, beginning right after the incident shock
    % Vs = array of specific volume, beginning right after the incident shock
    % P_as = pressure after incident shock
    % T_as = temperature after incident shock


function [ts, Vms, P_as, T_as]= particle_rre(mode, x, y, mach)
    %could be useful to introduce a "case" type of structure...
    if mode == 1
        H2_O2_conds;
    end
    
    %computation of pressure jump cooresponding to mach number selected
    xi = machToXi(mach,gamma_I,pi/2);
    P_as = P0*xi/101325; %conversion en atm
    
    %distance to point of incidence, when shock reaches particle
    r_0=y-x/tan(omega_rad); %has to be positive!
    if r_0<0
        "This code only follows particles in phase I!"
    else
        
        % should think about accelerating computation...
        %
        %
        
        %When in zone 1 between incident shock and expansion fan
        Msh=sqrt(xiToSqMach(xi,gamma_I,pi/2));
        Mi=Msh/sin(omega_rad);
        delta_1=atan(sqrt(tanDefSq(xi,Mi,gamma_I)));
        mu_1=asin(1/Mi); %velocity vector-expansion fan angle at beginning...
        %... of expansion fan
        t_1=1/sqrt(gamma_I*R_I*T0_I)*r_0*sin(omega_rad-delta_1+mu_1); %time spent...
        %between the incident shock and the expansion
        temp_jump=xiToTempJ(xi,gamma_I); %temperature jump of incident shock
        spec_vol_jump=((gamma_I-1)*xi+(gamma_I+1))/((gamma_I+1)*xi+(gamma_I-1));
        r_1=r_0/(Mi*sin(omega_rad-delta_1)); %radius when arriving at expansion


        %When in zone 2, behind the expansion fan
        Mr=sqrt(postShockMachSq(xi,Mi,gamma_I));%Mach when arriving at reflection

        %Computing expansion deflection and pressure jump
        Mt=sqrt((gamma_I*T0_I*mu_II)/(gamma_II*T0_II*mu_I))*Mi;
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
        T_post_exp=T0_I*temp_jump*exp_temp_jump; %temperature after expansion
        exp_spec_vol_jump=((1+(gamma_I-1)/2*Mr^2)...
            /(1+(gamma_I-1)/2*M_post_exp^2))^(1/(1-gamma_I)); %spec. vol. jump...
        %... after expansion
        Vm_post_exp=Vm0_I*spec_vol_jump*exp_spec_vol_jump;

        %In expansion fan
        no_dev_pts=1000; %number of points to solve differential equation
        ddelta=delta_2/(no_dev_pts); %deviation step

        T_as=T0_I*temp_jump; %temperature before expansion fan

        u1=Mr*sqrt(gamma_I*R_I*T_as); %speed behind expansion fan
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
        Texp=T_as; %temperature at current position
        Vmexp=Vm0_I*spec_vol_jump; %spec. vol. at current position

        for i=1:no_dev_pts
            dMsqexp=2*Msqexp*(1+(gamma_I-1)/2*Msqexp)/sqrt(Msqexp-1)*ddelta; %square ...
            %... Mach variation

            %Computing pressure and temperature
            dPexp=-Pexp*(gamma_I/2)/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %pressure variation
            Pexp=Pexp+dPexp; %new pressure
            dTexp=-Texp*((gamma_I-1)/2)/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %temperature variation
            Texp=Texp+dTexp; %new temperature
            dVmexp=Vmexp/2*1/(1+(gamma_I-1)/2*Msqexp)*dMsqexp; %spec. vol. variaton
            Vmexp=Vmexp+dVmexp;
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

        %construction of time array
        %Chemkin needs time to start after the incident shock
        
        refr_dur=experiment_dim/(Msh*sqrt(gamma_I*R_I*T0_I)); 

        % t_1 time spent between the incident shock and the expansion        
        t_after=tsexp(end)+refr_dur;
        
        n_s = round(t_1/3e-9);
        n_t = round(refr_dur/3e-9);
        ts=zeros(1,no_dev_pts+n_s+n_t);

        ts(1:n_s) = linspace(0,t_1,n_s);
        ts(n_s+1:n_s+no_dev_pts)=t_1+tsexp;
        ts(n_s+1+no_dev_pts:n_s+n_t+no_dev_pts)=linspace(ts(n_s+no_dev_pts),t_after,n_t);

        % construction of specific volume array
        Vms=zeros(1,no_dev_pts+n_s+n_t);
        Vms(1:n_s)= Vm0_I*spec_vol_jump;
        Vms(n_s+1:n_s+no_dev_pts)= Vmexps;
        Vms(n_s+no_dev_pts+1:n_s+no_dev_pts+n_t) = Vm_post_exp;
        
        Vms=Vms/Vm0_I;

        if ~(length(ts) == length(Vms))
            "PB DE DIMENSION"
        end
    end
end