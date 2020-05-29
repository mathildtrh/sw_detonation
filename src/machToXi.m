function xi=machToXi(Mach,gamma,phi)
    %calculates pressure jump when Mach is known
    %Inputs:
        %Mach: Mach number of oblique shock
        %gamma: ratio of constant heat
        % phi: shock - pre-shock flow angle 
    xi=((1-gamma)+2*gamma*(Mach*sin(phi))^2)/(1+gamma); %...
        %... transmitted pressure jump
end