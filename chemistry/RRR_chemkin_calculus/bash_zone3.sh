#!/bin/bash

cd $1
cp senk_zone3.inp ../senk.inp
cd ../
wine snkvtim.exe 
wine snkout.exe
cd -
cp ../t_T_ESP_Ther.dat ZN3.dat