% This MATLAB script is designed to get the Mach number of ignition
% for a RRE structure, depending on the particle studied

% main_RRE.m maybe needs to be run before the script, so that all the directory
% names are known and all the .inp files exist

% if main_RRE.m has already been ran, chemical calculus is not needed
% chemical_run has to be false
chemical_run = false;

%% Initialisation
disp("Running main_RRE.m ... It may take a while... ¯\_(ツ)_/¯");
run main_RRE.m;
disp("done");

ordinate = Y*100;
M_igni = zeros(1, length(Y));

%% Main loop
for part = 1:1:nP
    % mach is not reinitialized for each particle because we already
    % observed that it is a decreasing function of ordinate
    mach = 101;
    folder = dirnames(nP*(mach-1)+part);
    [~,temp]=plotChemkinRRE(folder);
    ignit = ignition(temp);
    % while the particle ignites, test a lower mach number
    while (ignit) && (mach>=1)
        mach=mach-1;
        folder = dirnames(nP*(mach-1)+part);
        [~, temp]=plotChemkinRRE(folder);
        ignit = ignition(temp);
    end
    mach = mach +1;
    if mach == 202
        M_igni(part) = mach_numbers(mach-1);
        disp(["Particle ", num2str(part), " seems to ignite before Mach = 7"])
    elseif mach == 1
        M_igni(part) = NaN;
        disp(["Particle ", num2str(part), " doesn't ignite before Mach = 9"])
    else
        M_igni(part) = mach_numbers(mach);
    end
end

M_flow = M_igni/sin(omega_deg*pi/180);

%% Plotting results
figure()
subplot(2,1,1)
plot(ordinate, M_flow,'*m')
hold on
xlabel("Ordinate (cm)")
ylabel("Mach number of ignition")
title ('Evolution of the Mach number of the flow which causes ignition');
grid on 
grid minor

subplot(2,1,2)
plot(ordinate, M_flow*sqrt(gamma_I*R_I*temp_I), '*r')
hold on
xlabel("Ordinate (cm)")
ylabel("Shock speed of ignition (m/s)")
title ('Evolution of the speed of the flow which causes ignition');
grid on 
grid minor

figure()
subplot(2,1,1)
plot(X*100, M_flow,'*m')
hold on
xlabel("Abscissa (cm)")
ylabel("Mach number of ignition")
title ('Evolution of the Mach number of the flow which causes ignition');
grid on 
grid minor

subplot(2,1,2)
plot(X*100, M_flow*sqrt(gamma_I*R_I*temp_I), '*r')
hold on
xlabel("Abscissa (cm)")
ylabel("Shock speed of ignition (m/s)")
title ('Evolution of the speed of the flow which causes ignition');
grid on 
grid minor