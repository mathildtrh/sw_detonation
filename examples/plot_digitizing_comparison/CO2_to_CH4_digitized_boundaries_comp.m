% Mathilde
% First program to compare plot digitizing
% based on graph found in Abd-El-Fattah's and L. F. Henderson 1978 paper
% Shock Waves at slow-fast gas interface, fig 13.

clear all

% Computing Yann's data (copied from
% examples/CO2_to_CH4/CO2_to_CH4_sys_lims.m)

fid=fopen('article_syst_lims.txt','r');
read_content=false;
no_graph_linesY=6;
points_Yann=zeros(2*no_graph_linesY,50);
line_index=0;
no_read_lines=0;
for graph_line_no=1:no_graph_linesY
    no_read_points=0;
    file_line=fgetl(fid);
    while ~contains(file_line,"Line")
        file_line=fgetl(fid);
    end
    if read_content
        no_read_lines=no_read_lines+1;
        line_index=line_index+1;
        file_line=fgetl(fid);
        while ~strcmp('',file_line) && ischar(file_line)
            point=str2num(file_line);
            no_read_points=no_read_points+1;
            points_Yann(2*no_read_lines-1:2*no_read_lines,no_read_points+1)=...
                point';
            file_line=fgetl(fid);
        end
        points_Yann(2*no_read_lines-1,1)=no_read_points;
    else
        read_content=true;
    end
end
fclose(fid);


% Computing Mathilde's data

fid=fopen('fig13-abdelfattah-henderson1978.txt','r');
%initialisation de la lecture du fichier
file_line=fgetl(fid);
while ~contains(file_line,"Line") %on va jusqu'à la 1e line qui contient "Line"
    file_line=fgetl(fid);
end
no_graph_linesM=5;
points_Math=zeros(2*no_graph_linesM,90);
line_index=0;
no_read_lines=0;
for graph_line_no=1:no_graph_linesM
    no_read_points=0;
    file_line=fgetl(fid); %cette ligne contient bien des chiffres? à vérifier + erreur
    %file_line
    no_read_lines=no_read_lines+1;
    line_index=line_index+1;
    file_line=fgetl(fid);
    while ~strcmp('',file_line) && ischar(file_line)
        [point, tf]=str2num(file_line);
        if ~tf
            print("Panique! Le point n'est pas un point")
        end 
        no_read_points=no_read_points+1;
        points_Math(2*no_read_lines-1:2*no_read_lines,no_read_points+1)=...
            point';
        file_line=fgetl(fid);
    end
    points_Math(2*no_read_lines-1,1)=no_read_points;
end
fclose(fid);

%Ploting limits + legend
%

cols=['m';'r';'g';'c';'b'];
for i=1:no_read_lines
    hold on
    no_read_points=points_Yann(2*i-1,1);
    plot(points_Yann(2*i-1,2:no_read_points+1),points_Yann(2*i,2:no_read_points+1),[cols(i) 'o']);
    %ploting Yann's digitizing
end

for j=1:no_read_lines
    hold on
    no_read_points=points_Math(2*j-1,1);
    plot(points_Math(2*j-1,2:no_read_points+1),points_Math(2*j,2:no_read_points+1),[cols(j) '*']);
    %ploting Mathilde's digitizing
end

legends={"RRE<->... Yann boundary points",...
    "RRR<->BPR Yann boundary points",...
    "BPR<->FNR Yann boundary points",...
    "FPR<->TNR Yann boundary points",...
    "TNR<->LSR Yann boundary points",...
    "RRE<->... Mathilde boundary points",...
    "RRR<->BPR Mathilde boundary points",...
    "BPR<->FNR Mathilde boundary points",...
    "FPR<->TNR Mathilde boundary points",...
    "TNR<->LSR Mathilde boundary points"};

lgd = legend(legends);
lgd.Location = "eastoutside";
xlabel("$\omega_i$ (deg)",'interpreter','latex')
ylabel("\chi")
title("Digitized boundaries for CO2->CH4 refraction, from fig 13, Adb-el-Fattah and Henderson 1978 paper")
xlim([25,90])
ylim([0,1])
hold off