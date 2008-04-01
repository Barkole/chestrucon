#!/bin/bash
echo "CONTCAR to INS"
echo "Name of Input-file?"
read INPUT
# INPUT=Test_CONTCAR
echo "Type in atom types in order as in CONTCAR/POTCAR!"
read ATOMS_1
echo $ATOMS_1 > ATOMS_1 
echo "Title for ins-file?"
read TITLE
echo "TITL " $TITLE > TITLE
cat ATOMS_1 $INPUT > AWK-INPUT
LC_ALL=C ./contcar2ins.awk AWK-INPUT > $INPUT.1
cat TITLE $INPUT.1 > $INPUT.ins
rm ATOMS_1 AWK-INPUT TITLE $INPUT.1
