% Make input for first CV calculation
nom='senk_zone2.inp';
fid = fopen(nom,'w');
fprintf(fid,'%s\n','CONV');
fprintf(fid,'%s','TIME', ' ');
fprintf(fid,'%i\n',t_2);
fprintf(fid,'%s','DELT', ' ');
fprintf(fid,'%i\n',1E-12);
fprintf(fid,'%s','REAC H2', ' ');
fprintf(fid,'%i\n',XH2);
fprintf(fid,'%s','REAC O2', ' ');
fprintf(fid,'%i\n',XO2);
fprintf(fid,'%s','PRES', ' ');
fprintf(fid,'%i\n',P2);
fprintf(fid,'%s','TEMP', ' ');
fprintf(fid,'%i\n',T2);
fprintf(fid,'%s','END');
fclose(fid);