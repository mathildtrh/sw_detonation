% Extract T, P, composition at the end of first CV calculation

[Header1, M1]=hfile('CV1.dat');
T=M1(end,2);
P=M1(end,3)/101325;
XOHs=M1(end,6);
XH=M1(end,7);
XH2=M1(end,8);
XO=M1(end,9);
XO2=M1(end,10);
XOH=M1(end,11);
XHO2=M1(end,12);
XH2O2=M1(end,13);
XH2O=M1(end,14);


% Make input for VTIM calculation
nom='senk_Exp.inp';
fid = fopen(nom,'w');
fprintf(fid,'%s\n','VTIM');
fprintf(fid,'%s','TIME', ' ');
fprintf(fid,'%i\n',tfinal2);
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