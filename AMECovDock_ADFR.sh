#!/bin/bash
# Given a .smi file with format: SMILE ID
# generates pdb files from that file
# generates pdbqt files from the created pdb files
# attempts to dock the pdbqt file to the provided receptor.trg file
#   chosen atoms are placed in output/ID/atoms.txt
# docked out.pdbqt files are placed in Res_Dock
# When the docking is finished, generates Res_Dock/scores.txt to display the docking scores
#
# USAGE: ./AMECovDock_ADFR.sh

# Programmer: Darin Hauner
###############  Modify These Values  #########################################
receptor="6wqf"
LIGANDS="Fin_test"
ADFRPATH="/people/jame118/ADFRsuite-1.0/ADFRsuite_x86_64Linux_1.0/bin/"
###############################################################################

PDBFILES='pdb/*.pdb'

##mkdir -m777 data
# to be automated
# agfr -r 6wqf.pdbqt -b user -17.00 -5.00 15.00 40.00 40.00 40.00 -c 1379 1382 -t 1377 -x A:CYS145 -o 6wqf_cov
rm -rf output
rm -rf Res_Dock

mkdir -m777 output
mkdir -m777 Res_Dock

# Run the reactions and generate the PDB files
python3 AMECovDock_Reactions.py

for f in $PDBFILES
do
    basepath=${f%.pdb}
    id=$(echo $basepath | awk -F'/' '{print $2}')
    base=$(echo $basepath | awk -F'_' '{print $(NF)}')
    echo "ID: " $id "BASE: " $base "BASEPATH: " $basepath
    mkdir output/${id}/

    # moves into and out of data to avoid file not found error when running prepare_ligand
    cd pdb/
    grep 'C   CYS A 145' ${id}.pdb > tmp
    grep 'CA  CYS A 145' ${id}.pdb >> tmp
    grep 'CB  CYS A 145' ${id}.pdb >> tmp
    grep 'SG  CYS A 145' ${id}.pdb >> tmp
    grep     'LIG A 999' ${id}.pdb >> tmp
    
    mv tmp ${id}.pdb

    ${ADFRPATH}prepare_ligand -R 1 -l ${id}.pdb
    mv ${id}.pdbqt ${receptor}_${id}.pdbqt

    ../find_group.sh ${receptor}_${id}.pdbqt
    group=$(cat found.txt)

    cd ..
    ${ADFRPATH}adfr -l pdb/${receptor}_${id}.pdbqt -t ${receptor}_cov.trg --jobName covalent -C $group --nbRuns 8 --maxEvals 100000 -O --seed 1

    mv ${receptor}_${id}_covalent_out.pdbqt Res_Dock/ 
    mv ${receptor}_${id}_covalent.dro ${receptor}_${id}_covalent_summary.dlg output/${id}/
    mv pdb/found.txt output/${id}/atoms.txt
done
#SCORE
for x in Res_Dock/*; do echo $x $(sed -n 3p $x) >> output.txt; done; sort -k4 -n output.txt > scores.txt
rm output.txt
mv scores.txt Res_Dock/
#DONE

