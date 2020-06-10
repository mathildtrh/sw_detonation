#!/bin/bash

cd $1
cp senk_Exp.inp ../senk.inp
cp VTIM.dat ../VTIM.dat
cd ../
wine snkvtim.exe 
wine snkout.exe
cd -
cp ../t_T_ESP_Ther.dat Exp.dat