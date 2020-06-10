#!/bin/bash

cd $1
cp senk_CV2.inp ../senk.inp
cd ../
wine snkvtim.exe 
wine snkout.exe
cd -
cp ../t_T_ESP_Ther.dat CV2.dat