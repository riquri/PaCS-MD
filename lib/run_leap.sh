#!/bin/bash

LIB_DIR=$(cd $(dirname $0); pwd)

cat << eof1 > ./get_center.in
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

mol = loadPDB initial.pdb

center mol

quit

eof1

tleap -f get_center.in

#The center is at: 33.85, 23.83, 25.42

IFS="," read -ra coordinates <<< $(grep "The center is at" leap.log | sed -E 's/.*The center is at: ([^,]+), ([^,]+), ([^,]+).*/\1,\2,\3/')

for i in 0 1 2 ; do
    coordinates[i]=$(echo "scale=2; -1*${coordinates[i]}" | bc)
done

cat << eof1 > ./leap.in
#Use Amber ff14SB force field
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

list

mol = loadPDB initial.pdb

center mol

set mol box {50 50 50 }
translate mol { ${coordinates[0]} ${coordinates[1]} ${coordinates[2]} }
solvateBox mol TIP3PBOX 0.0

addIons2 mol Cl- 0
addIons2 mol Na+ 0
charge mol

saveAmberParm mol ./leap.prmtop ./leap.inpcrd

savePDB mol leap.pdb

quit

eof1

tleap -f leap.in

$LIB_DIR/run_acpype.sh



