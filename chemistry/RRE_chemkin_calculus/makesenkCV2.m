% Extract T, P, composition at the end of first CV calculation

[Header2, M2]=hfile('Exp.dat');
T=M2(end,2);
P=M2(end,3)/101325;
XOHs=M2(end,6);
XH=M2(end,7);
XH2=M2(end,8);
XO=M2(end,9);
XO2=M2(end,10);
XOH=M2(end,11);
XHO2=M2(end,12);
XH2O2=M2(end,13);
XH2O=M2(end,14);


% Make input for second CV calculation
nom='senk_CV2.inp';
fid = fopen(nom,'w');
fprintf(fid,'%s\n','CONV');
fprintf(fid,'%s','TIME', ' ');
fprintf(fid,'%i\n',tfinal3);
fprintf(fid,'%s','DELT', ' ');
fprintf(fid,'%i\n',1E-12);
fprintf(fid,'%s','REAC OH*', ' ');
fprintf(fid,'%i\n',XOHs);
fprintf(fid,'%s','REAC H', ' ');
fprintf(fid,'%i\n',XH);
fprintf(fid,'%s','REAC H2', ' ');
fprintf(fid,'%i\n',XH2);
fprintf(fid,'%s','REAC O', ' ');
fprintf(fid,'%i\n',XO);
fprintf(fid,'%s','REAC O2', ' ');
fprintf(fid,'%i\n',XO2);
fprintf(fid,'%s','REAC OH', ' ');
fprintf(fid,'%i\n',XOH);
fprintf(fid,'%s','REAC HO2', ' ');
fprintf(fid,'%i\n',XHO2);
fprintf(fid,'%s','REAC H2O2', ' ');
fprintf(fid,'%i\n',XH2O2);
fprintf(fid,'%s','REAC H2O', ' ');
fprintf(fid,'%i\n',XH2O);
fprintf(fid,'%s','PRES', ' ');
fprintf(fid,'%i\n',P);
fprintf(fid,'%s','TEMP', ' ');
fprintf(fid,'%i\n',T);
fprintf(fid,'%s','END');
fclose(fid);