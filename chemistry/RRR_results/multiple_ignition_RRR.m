% This MATLAB script is designed to compare the several results obtained
% for a RRR structure, depending on the mach number selected or the
% particle chosen

% main_RRR.m maybe needs to be run before the script, so that all the directory
% names are known and all the .inp files exist

% if main_RRR.m has already been ran, chemical calculus is not needed
% chemical_run has to be false
chemical_run = false;

disp("Running main_RRR.m ... It may take a while... ¯\_(ツ)_/¯");
run main_RRR.m;
disp("done");
mode = 4;

%% mode = 1 : one particle // all mach numbers
if mode == 1
    part = 2; % CHANGE
    for mach = 1:2:nM-1
        folder = dirnames(nP*(mach-1)+part);
        lgd = strcat("Mach = ", num2str(mach_flow(mach)));
        [time,temp]=plotChemkinRRR(folder);
        p = plot(time, temp, 'DisplayName', lgd);
        p.LineWidth = 1.2;
        legend('Location', 'eastoutside')
        hold on
        xlabel("Time (s)")
        ylabel("Temperature (K)")
        title ({"Evolution of temperature with time,"; ...
            "for several incident shock strengths,"; ...
            strcat("fluid particle : x = ", num2str(X(part)*100), "cm")});
    end
end

%% mode = 2 : one figure for each particle
if mode == 2
    for part = 1:1:nP
        figure()
        for mach = 1:1:nM
            folder = dirnames(nP*(mach-1)+part);
            lgd = strcat("Mach = ", num2str(mach_flow(mach)));
            [time,temp]=plotChemkinRRR(folder);
            p = plot(time, temp, 'DisplayName', lgd);
            p.LineWidth = 1.2;
            legend('Location', 'eastoutside')
            hold on
            xlabel("Time (s)")
            ylabel("Temperature (K)")
            title ({"Evolution of temperature with time,"; ...
                "for several incident shock strengths,"; ...
                strcat("fluid particle : x = ", num2str(X(part)*100), "cm")});
        end
    end
end

%% mode = 3 : one mach number // all particles
if mode == 3
    mach = 12; % CHANGE
    for part = 1:1:nP
        folder = dirnames(nP*(mach-1)+part);
        lgd = strcat("X = ", num2str(X(part)));
        [time,temp]=plotChemkinRRR(folder);
        p = plot(time, temp, 'DisplayName', lgd);
        p.LineWidth = 1.2;
        legend('Location', 'eastoutside')
        hold on
        xlabel("Time (s)")
        ylabel("Temperature (K)")
        title ({"Evolution of temperature with time,"; ...
            "for several fluid particles,";...
            strcat("Mach number of shock = ", num2str(mach_flow(mach)))});
    end
end

%% mode = 4 : one figure for each mach number
if mode == 4
    for mach = [10 20 30]
        figure()
        for part = 1:1:nP
            folder = dirnames(nP*(mach-1)+part);
            lgd = strcat("X = ", num2str(X(part)));
            [time,temp]=plotChemkinRRR(folder);
            p = plot(time, temp, 'DisplayName', lgd);
            p.LineWidth = 1.2;
            legend('Location', 'eastoutside')
            hold on
            xlabel("Time (s)")
            ylabel("Temperature (K)")
            title ({"Evolution of temperature with time,"; ...
                "for several fluid particles,"; ...
                strcat("Mach number of shock = ", num2str(mach_flow(mach)))});
        end
    end
end