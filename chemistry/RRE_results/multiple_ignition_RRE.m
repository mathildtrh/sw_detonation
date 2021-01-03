% This MATLAB script is designed to compare the several results obtained
% for a RRE structure, depending on the mach number selected or the
% particle chosen

% main_RRE.m maybe needs to be run before the script, so that all the directory
% names are known and all the .inp files exist

% if main_RRE.m has already been ran, chemical calculus is not needed
% chemical_run has to be false
chemical_run = false;

disp("Running main_RRE.m ... It may take a while... ¯\_(ツ)_/¯");
run main_RRE.m;
disp("done");

cd '/home/mathilde/Documents/PRe_Detonation/sw_detonation/chemistry/RRE_chemkin_calculus';
mode = 2;

%% mode = 1 : one particle // all mach numbers
if mode == 1
    part = 1; % CHANGE
    for mach = 1:1:nM
        folder = dirnames(nP*(mach-1)+part);
        lgd = strcat("Mach = ", num2str(round(mach_flow(mach),1)));
        [time,temp]=plotChemkinRRE(folder);
        p = plot(time, temp, 'DisplayName', lgd);
        p.LineWidth = 2;
        legend('Location', 'east')
        hold on
        xlabel("Time (s)",'FontSize', 18)
        ylabel("Temperature (K)", 'FontSize', 18)
        %title ({"Evolution of temperature with time,"; ...
            %"for several incident shock strengths,"; ...
            %strcat("fluid particle : x = ", num2str(X(part)*100), "cm")});
    end
end

%% mode = 2 : one figure for each particle
if mode == 2
    tosave = zeros(50000, 2*nP*nM);
    for part = 1:1:nP
        %figure()
        for mach = 1:1:nM
            folder = dirnames(nP*(mach-1)+part);
            lgd = strcat("Mach = ", num2str(mach_flow(mach)));
            [time,temp]=plotChemkinRRE(folder);
            len = length(time);
            %p = plot(time, temp, 'DisplayName', lgd);
            tosave(1:len,nP*(mach-1)+part)=time';
            tosave(1:len,nP*nM+nP*(mach-1)+part)=temp';
            %legend('Location', 'east')
            %hold on
            %xlabel("Time (s)")
            %ylabel("Temperature (K)")
            %title ({"Evolution of temperature with time,"; ...
                %"for several incident shock strengths,"; ...
                %strcat("fluid particle : x = ", num2str(X(part)*100), "cm")});
        end
    end
    readme = "Particule 1, mach 1 : time en colonne 1, temp en colonne nP*nM +1";
    cd('/home/mathilde/Documents/PRe_Detonation/sw_detonation/chemistry/RRE_chemkin_calculus');
    save('ignition_RRE');
end

%% mode = 3 : one mach number // all particles
if mode == 3
    mach = 1; % CHANGE
    for part = 1:1:nP
        folder = dirnames(nP*(mach-1)+part);
        lgd = strcat("X = ", num2str(X(part)));
        [time,temp]=plotChemkinRRE(folder);
        p = plot(time, temp, 'DisplayName', lgd);
        legend('Location', 'east')
        hold on
        xlabel("Time (s)")
        ylabel("Temperature (K)")
        %title ({"Evolution of temperature with time,"; ...
            %"for several fluid particles,";...
            %strcat("Mach number of shock = ", num2str(mach_flow(mach)))});
    end
end

%% mode = 4 : one figure for each mach number
if mode == 4
    for mach = [12 23 25]
        figure()
        for part = 1:1:nP
            folder = dirnames(nP*(mach-1)+part);
            lgd = strcat("X = ", num2str(X(part)*100), 'cm');
            [time,temp]=plotChemkinRRE(folder);
            p = plot(time, temp, 'DisplayName', lgd);
            legend('Location', 'east')
            p.LineWidth = 1.3;
            hold on
            grid on
            xlabel("Time (s)")
            ylabel("Temperature (K)")
            title ({"Evolution of temperature with time,"; ...
                "for several fluid particles,"; ...
                strcat("Mach number of shock = ", num2str(mach_flow(mach)))});
        end
    end
end