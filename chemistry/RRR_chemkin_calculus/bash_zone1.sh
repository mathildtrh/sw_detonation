#!/bin/bash

cd $1
cp senk_zone1.inp ../senk.inp
cd ../
wine snkvtim.exe 
wine snkout.exe
cd -
cp ../t_T_ESP_Ther.dat ZN1.dat