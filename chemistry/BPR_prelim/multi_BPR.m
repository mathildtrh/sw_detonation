% This Matlab script is designed to explore a rectangle in the omega-chi
% plane to determine whether several BPR forms can coexist in their
% definition domain

% first BPR with reflected shock
% second BPR with reflected expansion

mode = 2;
% mode 1 : H2-O2 // He system
% mode 2 : CO2 // CH4 system
if mode == 1
    load H2_O2_to_He_limitsMathilde.mat
    run H2_O2_conds.m
    omegas = linspace(20,60);
else
    load CO2_to_CH4_sys_limits_data.mat
    gamma_I=1.288; 
    gamma_II=1.303; 
    mu_I=44.01;
    mu_II=16.04;
    R_u=8.314;
    R_I=R_u/(mu_I*10^-3);
    R_II=R_u/(mu_II*10^-3);
    omegas = linspace(30,80);
end

cst = [gamma_I, gamma_II, R_I, R_II];

len_om = length(omegas);
chis = linspace(0,1);
len_ch = length(chis);
result = -1*ones(len_om, len_ch);

% result(om,ch) = 1 <=> at the point omegas(om);chis(ch), BPR has a shock
% result(om,ch) = 0 <=> at the point omegas(om);chis(ch), BPR has an expansion
% result(om,ch) = -1 <=> at the point omegas(om);chis(ch), it is not BPR

for om = 1:1:len_om
    for ch = 1:1:len_ch
        % boolean to know whether calculation is needed
        calc = true;
        if chis(ch)>chis_RRR_BPR(1)
            % is it RRR?
            [~, index]=min(abs(chis_RRR_BPR-chis(ch)));
            if omegas(om)<omegas_RRR_BPR(index)
                calc=false;
            end
            
        else
            % is it RRE?
            [~, index]=min(abs(chis_RRE-chis(ch)));
            if omegas(om)<omegas_RRE(index)
                calc=false;
            end
        end
        % is it FNR?
        [~, index]=min(abs(chis_BPR_FNR-chis(ch)));
        if omegas(om)>omegas_BPR_FNR(index)
            calc=false;
        end
        
        % If we are in the BPR domain, then calculation is needed
        if calc
            mach = sqrt(xiToSqMach(1/chis(ch), gamma_I,pi/2));
            result(om, ch) = BPRdisc(cst, mach, omegas(om), 101325, 600, 1138);
        end
    end
end

mycolormap = [1 1 1; 0.9 0.9 0 ; 0 0.9 0.9];
figure()
hold on
pcolor(omegas, chis, result')
xlabel('\omega (deg)')
ylabel('\chi')
title('Map of configurations for BPR structure')
colormap(mycolormap);