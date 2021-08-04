#/bin/bash/
# finds the three anchor atoms in a pdbqt file and writes them to found.txt
# programmer: Darin Hauner
val1=$(grep 'C   CYS' $1 | grep 'ATOM' | awk '{print $2}')
val2=$(grep 'CA' $1 | grep 'ATOM' | awk '{print $2}')
val3=$(grep 'CB' $1 | grep 'ATOM' | awk '{print $2}')
echo $val1 $val2 $val3 > found.txt
