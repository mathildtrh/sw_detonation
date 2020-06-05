% This program computes the evolution of normalized specific volume
% of several lagrangian particles
% for different shock mach and a fixed angle of incidence

clear all
close all

mode = 1;

mach_numbers = [4];

% X = [-0.9, 0.5, 1, 1.5, -0.5];
% Y = [3, 3, 7, 11, 12];
X = [0.5, 1, 1.5];
Y = [3, 5, 8];
nM = length(mach_numbers);
nP = length(X);

%filenames = [""]*nM*nP;
% names repartition
% mach(1)_x1_y1, mach(1)_x2_y2, mach(1)_x3_y3,...
% ... mach(2)_x1_y1, mach(2)_x2_y2, mach(2)_x3_y3,...

for m = 1:nM
    for p = 1:nP
        filenames(nP*(m-1)+p) = strcat("chemistry/vtim_matlab/mach", num2str(mach_numbers(m)),...
            "_x", num2str(X(p)), "_y", num2str(Y(p)));
        senknames(nP*(m-1)+p) = strcat("chemistry/senk_matlab/mach", num2str(mach_numbers(m)),...
            "_x", num2str(X(p)), "_y", num2str(Y(p)));
    end
end

X=X*1e-2;
Y=Y*1e-2;

for j=1:nM
    for i = 1:nP
        [t, V, P1, T1] = particle_rre(mode, X(i), Y(i),mach_numbers(j));
        "Done"
        file = fopen(filenames(nP*(j-1)+i),'w');
        fprintf(file, '%1.4E   %1.5f\n', cat(1,t(1:end-1),V(1:end-1)));
        fprintf(file, '%1.4E   %1.5f', t(end),V(end));
        fclose(file);
        file2 = fopen(senknames(nP*(j-1)+i),'w');
        fprintf(file2, '%s\n', 'VTIM');
        fprintf(file2, '%s %E\n', 'TIME', t(end));
        fprintf(file2, '%s %e\n', 'DELT', 1e-12);
        fprintf(file2, '%s %s  %1.4f\n', 'REAC', 'H2', 2/3);
        fprintf(file2, '%s %s  %1.4f\n', 'REAC', 'O2', 1/3);
        fprintf(file2, '%s %2.2f\n', 'PRES', round(P1,2));
        fprintf(file2, '%s %4.0f\n', 'TEMP', round(T1));
        fprintf(file2, '%s', 'END');
        fclose(file2);
    end
end