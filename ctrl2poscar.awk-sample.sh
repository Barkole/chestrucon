#!/bin/bash
echo "Convert Ca3N2.CTRL to POSCAR"
LC_ALL=C ./ctrl2poscar.awk Ca3N2.CTRL > Ca3N2.POSCAR
