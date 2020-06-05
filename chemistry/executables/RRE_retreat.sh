#!/bin/bash

# Script d'automatisation du traitement des données MATLAB
# pour l'évolution d'une particule lagragienne
# Initialement développé pour les structures RRE

echo "Bienvenue dans le script de calculs chimiques !"


# Vérifier qu'on est dans le bon dossier
# celui qui contient les .exe

#check=`ls *.exe`
#if [ -z $check ]; then
#    echo "Pas de fichiers executable"
#    break
#fi

# Récupérer tous les noms de fichier intéressants
# dans le bon répertoire, passé en argument 1
list_VTIM=`ls $1`

# faire tourner l'interpreteur 
echo -e "Execution de l'interpreteur CHEMKIN... \n"
wine INTERP.exe
# vérifier execution correcte dans chem.out
tail -6 chem.out | head -1
echo -e "\n\n"

# Pour chaque fichier intéressant, faire
for file in $list_VTIM; do
    echo "-----------> " $file "  --- Traitement en cours --------- "
    # copier le fichier VTIM dans VTIM.dat
    cp ../vtim_matlab/$file VTIM.dat
    cp ../senk_matlab/$file senk.inp

    # exécuter snkvtim avec wine
    wine snkvtim.exe
    echo "snkvtim done"

    # exécuter snkout avec wine
    wine snkout.exe
    echo "snkout done"
    # copier t_T_ESP_Ther.dat dans le bon fichier
    cp t_T_ESP_Ther.dat ../output_chemkin/
    mv ../output_chemkin/t_T_ESP_Ther.dat ../output_chemkin/$file.dat
    echo "copy = ok"
done

#for i in ../output_chemkin/VTIM_*;do 
#mv $i ../output_chemkin/out_${i#../output_chemkin/VTIM_}
#done
echo "Fin du programme !"