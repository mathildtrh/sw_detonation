%% WARNING : this function has been changed a certain amount of times
% it needs to be read carefully because some lines of code are to be
% commented/uncommented depending on the way you want the code to run

function [before_transition,Msj,Mst,Msi]=fromFPRToTNR(xi_i,omega_rad,gamma_I,...
    gamma_II,mu_I,mu_II,ny,varargin)
    %Determines if FPR to TNR transition was made
    %Inputs:
        %xi_i: incident shock wave pressure jump
        %omega_rad: angle of incidence (rad)
        %gamma_I, gamma_II: ratios of constant heat
        %mus: molecular masses
        %ny: xi-axis resolution when getting polars
        %vargarin{1}, solve_Msj_eq: bool determining if Msj 
            %... equation needs to be solved or not, used to reduce...
            %... computation time
        %varargin{2}, Msj: Msj value, used if solve_Msj_eq==false
        %varargin{3}, temp_ratio: ratio of temperatures=T_I/T_II
    clear a
    solve_Msj_eq=true;
    temp_ratio=1;
    if nargin>=8
        solve_Msj_eq=varargin{1};
        if nargin>=10
            temp_ratio=varargin{3};
        end
        if nargin>=11
            x1=varargin{4};
            x2=varargin{5};
            x3=varargin{6};
            x4=varargin{7};
            x5=varargin{8};
        end
    end
    
    a_I = sqrt(gamma_I/mu_I*temp_ratio);
    a_II = sqrt(gamma_II/mu_II);
    
    
    Msi=sqrt(xiToSqMach(xi_i,gamma_I,pi/2)); %incident shock Mach
    Mi=Msi/sin(omega_rad); %incident free stream Mach
    Vi = Msi*a_I;

    %Computing Mst:
    
%     %1st correction (constant coeff)
%     x1 = 0.98;
%     b=(gamma_II+1)/(gamma_I+1)*(Vi^2-a_I^2)/Vi * x1;
    
%     %polynomial correction
%     x1 = -0.3974;
%     x2 = -0.0369;
%     x3 = 0.8176;
%     b=(gamma_II+1)/(gamma_I+1)*(Vi^2-a_I^2)/Vi *(x1*omega_rad^2 +x2*omega_rad +x3);
    
%     %trigo correction
%     x1 = 0.7539;
%     b=(gamma_II+1)/(gamma_I+1)*(Vi^2-a_I^2)/Vi *(x1*cos(omega_rad));

    %multitrigo correction
    b=(gamma_II+1)/(gamma_I+1)*(Vi^2-a_I^2)/Vi *(x1*cos(x4*omega_rad)+x2*sin(x5*omega_rad)+x3);
    
    Vt=1/2*(b+sqrt(b^2+4*a_II^2)); %second and third term are for correction
    Mst = Vt/a_II;
    xi_t=((1-gamma_II)+2*gamma_II*Mst^2)/(1+gamma_II); %transmitted pressure jump
    
    Mt=Mi*a_I/a_II;
    phi_t = asin(Mst/Mt);
    t_dev=postShockDeflection(Mt,gamma_II,phi_t); % post shock deflection j
    
    %computing j precursor shock Mach
    if solve_Msj_eq
        syms Msj
        eqn= ((1-gamma_I)+2*gamma_I*Msj^2)/(1+gamma_I)...
            ==xi_t*(1-(gamma_II-1)/(gamma_I+1)...
                *sqrt((gamma_I*mu_II)/(gamma_II*mu_I))*...
            (Msj-1/Msj))^(2*gamma_II/(gamma_II-1));
        a=[vpasolve(eqn,Msj,[1,Inf])];
    else
        a=[varargin{2}];
    end
    if isempty(a) %if the equation has no solution : Msj<1 --> still FPR
        before_transition=true;
        Msj = 0.5;
    else

        phi_i = asin(Msi/Mi);
        i_dev=postShockDeflection(Mi,gamma_I,phi_i); % deflection of flow behind incident shock
        
        M1r=sqrt(postShockMachSq(xi_i,Mi,gamma_I)); % Mach number behind incident shock
        min_r_dev=i_dev-deltaMax(M1r,gamma_I); %minimum of r deviation
        
        if exist('a')==0
            disp("a doesn't exist")
        end
        Msj=a(1);
        if ~isnumeric(Msj)
            Msj=double(Msj);
        end
        
        xi_j=((1-gamma_I)+2*gamma_I*Msj^2)/(1+gamma_I); % pressure jump through j shock
        
        Mjn = sqrt(xiToSqMach(xi_j,gamma_I,pi/2)); % pre-shock mach j
        phi_j = asin(Mjn/Mi);
        j_dev=postShockDeflection(Mi,gamma_I,phi_j); % post shock deflection j
        
        M1k=sqrt(postShockMachSq(xi_j,Mi,gamma_I)); % Mach number behind j shock
        max_k_dev=-j_dev+deltaMax(M1k,gamma_I); %max of k deviation
        if max_k_dev<min_r_dev
            before_transition=false;
        else
            % plot to have some graphical help
            %[xis_i, deltas_i]=getPolarMathilde(Mi, gamma_I, ny, 0);
            %[xis_t, deltas_t]=getPolarMathilde(Mst, gamma_II, ny, 0);
            
            %[xis_j,deltas_j]=getPolarMathilde(Msj,gamma_I,ny,2,xi_t,t_dev);
            [xis_k,deltas_k]=getPolarMathilde(M1k,gamma_I,ny,2,xi_j,j_dev);
            [xis_r,deltas_r]=getPolarMathilde(M1r,gamma_I,ny,2,xi_i,i_dev);
            
%             figure()
%             hold on
%             semilogy(deltas_i*180/pi, xis_i, 'r', 'DisplayName', 'Incident')
%             semilogy(deltas_t*180/pi, xis_t, 'g', 'DisplayName', 'Transmitted')
%             semilogy(deltas_j*180/pi, xis_j, 'y', 'DisplayName', 'j shock')
%             semilogy(deltas_k*180/pi, xis_k, 'c', 'DisplayName', 'k shock')
%             semilogy(deltas_r*180/pi, xis_r, 'm', 'DisplayName', 'r shock')
%             legend();
%             grid on
%             grid minor
            
            i=ny+2;
            xi_k=xis_k(i); delta_k=deltas_k(i);
            [xi_coors,delta_coors,inds]=getPolarPoint(xis_r,deltas_r,0,xi_k);
            xi_r=xi_coors(1); delta_r=delta_coors(1);
            while (i<=2*ny+1) && ((xi_r==-1) || (delta_k<delta_r))...
                    && (xi_k>=xi_i)
                i=i+1;
                xi_k=xis_k(i);
                delta_k=deltas_k(i);
                [xi_coors,delta_coors,inds]=...
                    getPolarPoint(xis_r,deltas_r,0,xi_k);
                xi_r=xi_coors(1);delta_r=delta_coors(1);
            end
            if (i==2*ny+1) || (xi_k<xi_i)
                before_transition=false;
            else
                before_transition=true;
            end
        end
    end
end