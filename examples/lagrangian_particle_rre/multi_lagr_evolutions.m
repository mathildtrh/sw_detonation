% Program designed to plt several lagrangian particles evolution at the
% same time
clear all
close all

mode = 2;
X = [0.5, 1, 1.5]*1e-2;
Y = [3, 5, 8]*1e-2;
colors = ["b","r","g"];

blue = strcat("X = ", num2str(X(1)), ", Y = ", num2str(Y(1)));
red = strcat("X = ", num2str(X(2)), ", Y = ", num2str(Y(2)));
green = strcat("X = ", num2str(X(3)), ", Y = ", num2str(Y(3)));

legends = {blue, red, green};

for i = 1:3
    [t, P, T, V, ttl] = particle_evol(mode, X(i), Y(i));
    "Done"
    
    subplot(3,1,1)
    plot(t,P,colors(i)) %plotting pressure of particle
    title({ttl; 'Pressure evolution'})
    xlabel("time (micros)")
    ylabel("Pressure (bar)")
    legend(legends,'Location','southeast')
    
    hold on
    
    subplot(3,1,2)
    plot(t,T,colors(i)) %plotting temperature of particle
    title('Temp evolution (K)')
    xlabel("time (micros)")
    ylabel("Temperature (K)")
    legend(legends,'Location','southeast')
    
    hold on
    
    subplot(3,1,3)
    plot(t,V,colors(i)) %plotting temperature of particle
    title('Spec. vol. evolution')
    xlabel("time (micros)")
    ylabel("Normalized specific volume")
    ylim([0,1.1])
    legend(legends,'Location','southeast')
    
    hold on
end
