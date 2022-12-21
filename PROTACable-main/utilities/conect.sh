#!/bin/bash

for poi in poi*rec.pdb;

do
v1=$(echo "$poi" | cut -d '_' -f 2)
v2=$(echo "$poi" | cut -d '_' -f 3)
v3=${v1}_${v2}

for e3 in e3*_${v3}_*pdb;
do

e3=$e3
done
linker="${v3}_docked_short.pdb"
con="${v3}_cnct.idx"

## amending pdb for linker ligation ## LINKER FIRST ##
line=$(head -n 1 $con)
C1=$(cut -d',' -f1 <<<$line)
C2=$(cut -d',' -f2 <<<$line)
C3=$(cut -d',' -f3 <<<$line)
C4=$(cut -d',' -f4 <<<$line)
C5=$(cut -d',' -f5 <<<$line)
C6=$(cut -d',' -f6 <<<$line)

echo "${poi},${e3},${linker},${con},${C1},${C2},${C3},${C4},${C5},${C6}"
# poi starting point is 90
# first linker attachment is 91
# second linker attachment is 92
# third linker attachment is 93 (double bond)
# termination of a linker/no linker is 94
# e3 connection is 95

if [ -z "$C2" ];
then
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM $(($C1+1))  N   UNL/HETATM $(($C1+1))  N94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM  $(($C1+1))  N   UNL/HETATM  $(($C1+1))  N94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM   $(($C1+1))  N   UNL/HETATM   $(($C1+1))  N94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM    $(($C1+1))  N   UNL/HETATM    $(($C1+1))  N94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM $(($C1+1))  C   UNL/HETATM $(($C1+1))  C94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM  $(($C1+1))  C   UNL/HETATM  $(($C1+1))  C94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM   $(($C1+1))  C   UNL/HETATM   $(($C1+1))  C94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM    $(($C1+1))  C   UNL/HETATM    $(($C1+1))  C94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM $(($C1+1))  O   UNL/HETATM $(($C1+1))  O94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM  $(($C1+1))  O   UNL/HETATM  $(($C1+1))  O94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM   $(($C1+1))  O   UNL/HETATM   $(($C1+1))  O94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM    $(($C1+1))  O   UNL/HETATM    $(($C1+1))  O94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM $(($C1+1))  P   UNL/HETATM $(($C1+1))  P94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM  $(($C1+1))  P   UNL/HETATM  $(($C1+1))  P94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM   $(($C1+1))  P   UNL/HETATM   $(($C1+1))  P94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM    $(($C1+1))  P   UNL/HETATM    $(($C1+1))  P94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM $(($C1+1))  S   UNL/HETATM $(($C1+1))  S94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM  $(($C1+1))  S   UNL/HETATM  $(($C1+1))  S94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM   $(($C1+1))  S   UNL/HETATM   $(($C1+1))  S94 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM    $(($C1+1))  S   UNL/HETATM    $(($C1+1))  S94 UNL/" $poi
line1=$(awk "/UNL/{c++} c==$(($C1+1)){print NR;exit}" $poi)
C9=$(awk "NR==$line1 { print $3 }" $poi)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N94 /" $poi
sed -i "${line1}s/  N${init} /  N94 /" $poi
sed -i "${line1}s/  C${init}  /  C94 /" $poi
sed -i "${line1}s/  C${init} /  C94 /" $poi
sed -i "${line1}s/  O${init}  /  O94 /" $poi
sed -i "${line1}s/  O${init} /  O94 /" $poi
sed -i "${line1}s/  P${init}  /  P94 /" $poi
sed -i "${line1}s/  P${init} /  P94 /" $poi
sed -i "${line1}s/  S${init}  /  S94 /" $poi
sed -i "${line1}s/  S${init} /  S94 /" $poi

line1=$(awk "/UNL/{c++} c==$(($C1+1)){print NR;exit}" $poi)
C9=$(awk "NR==$line1 { print $2 }" $poi)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N94 /" $poi
sed -i "${line1}s/  N${init} /  N94 /" $poi
sed -i "${line1}s/  C${init}  /  C94 /" $poi
sed -i "${line1}s/  C${init} /  C94 /" $poi
sed -i "${line1}s/  O${init}  /  O94 /" $poi
sed -i "${line1}s/  O${init} /  O94 /" $poi
sed -i "${line1}s/  P${init}  /  P94 /" $poi
sed -i "${line1}s/  P${init} /  P94 /" $poi
sed -i "${line1}s/  S${init}  /  S94 /" $poi
sed -i "${line1}s/  S${init} /  S94 /" $poi
else
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM $(($C1+1))  N   UNL/HETATM $(($C1+1))  N90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM  $(($C1+1))  N   UNL/HETATM  $(($C1+1))  N90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM   $(($C1+1))  N   UNL/HETATM   $(($C1+1))  N90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM    $(($C1+1))  N   UNL/HETATM    $(($C1+1))  N90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM $(($C1+1))  C   UNL/HETATM $(($C1+1))  C90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM  $(($C1+1))  C   UNL/HETATM  $(($C1+1))  C90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM   $(($C1+1))  C   UNL/HETATM   $(($C1+1))  C90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM    $(($C1+1))  C   UNL/HETATM    $(($C1+1))  C90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM $(($C1+1))  O   UNL/HETATM $(($C1+1))  O90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM  $(($C1+1))  O   UNL/HETATM  $(($C1+1))  O90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM   $(($C1+1))  O   UNL/HETATM   $(($C1+1))  O90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM    $(($C1+1))  O   UNL/HETATM    $(($C1+1))  O90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM $(($C1+1))  P   UNL/HETATM $(($C1+1))  P90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM  $(($C1+1))  P   UNL/HETATM  $(($C1+1))  P90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM   $(($C1+1))  P   UNL/HETATM   $(($C1+1))  P90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM    $(($C1+1))  P   UNL/HETATM    $(($C1+1))  P90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM $(($C1+1))  S   UNL/HETATM $(($C1+1))  S90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM  $(($C1+1))  S   UNL/HETATM  $(($C1+1))  S90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM   $(($C1+1))  S   UNL/HETATM   $(($C1+1))  S90 UNL/" $poi
sed -i -z -e "s/UNL/UNL/$(($C1+1))" -e "s/HETATM    $(($C1+1))  S   UNL/HETATM    $(($C1+1))  S90 UNL/" $poi
line1=$(awk "/UNL/{c++} c==$(($C1+1)){print NR;exit}" $poi)
C9=$(awk "NR==$line1 { print $3 }" $poi)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N90 /" $poi
sed -i "${line1}s/  N${init} /  N90 /" $poi
sed -i "${line1}s/  C${init}  /  C90 /" $poi
sed -i "${line1}s/  C${init} /  C90 /" $poi
sed -i "${line1}s/  O${init}  /  O90 /" $poi
sed -i "${line1}s/  O${init} /  O90 /" $poi
sed -i "${line1}s/  P${init}  /  P90 /" $poi
sed -i "${line1}s/  P${init} /  P90 /" $poi
sed -i "${line1}s/  S${init}  /  S90 /" $poi
sed -i "${line1}s/  S${init} /  S90 /" $poi

line1=$(awk "/UNL/{c++} c==$(($C1+1)){print NR;exit}" $poi)
C9=$(awk "NR==$line1 { print $2 }" $poi)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N90 /" $poi
sed -i "${line1}s/  N${init} /  N90 /" $poi
sed -i "${line1}s/  C${init}  /  C90 /" $poi
sed -i "${line1}s/  C${init} /  C90 /" $poi
sed -i "${line1}s/  O${init}  /  O90 /" $poi
sed -i "${line1}s/  O${init} /  O90 /" $poi
sed -i "${line1}s/  P${init}  /  P90 /" $poi
sed -i "${line1}s/  P${init} /  P90 /" $poi
sed -i "${line1}s/  S${init}  /  S90 /" $poi
sed -i "${line1}s/  S${init} /  S90 /" $poi
fi


if [ -z "$C2" ];
then
:
else
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM $(($C2+1))  N   UNL/HETATM $(($C2+1))  N91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM  $(($C2+1))  N   UNL/HETATM  $(($C2+1))  N91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM   $(($C2+1))  N   UNL/HETATM   $(($C2+1))  N91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM    $(($C2+1))  N   UNL/HETATM    $(($C2+1))  N91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM $(($C2+1))  C   UNL/HETATM $(($C2+1))  C91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM  $(($C2+1))  C   UNL/HETATM  $(($C2+1))  C91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM   $(($C2+1))  C   UNL/HETATM   $(($C2+1))  C91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM    $(($C2+1))  C   UNL/HETATM    $(($C2+1))  C91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM $(($C2+1))  O   UNL/HETATM $(($C2+1))  O91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM  $(($C2+1))  O   UNL/HETATM  $(($C2+1))  O91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM   $(($C2+1))  O   UNL/HETATM   $(($C2+1))  O91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM    $(($C2+1))  O   UNL/HETATM    $(($C2+1))  O91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM $(($C2+1))  P   UNL/HETATM $(($C2+1))  P91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM  $(($C2+1))  P   UNL/HETATM  $(($C2+1))  P91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM   $(($C2+1))  P   UNL/HETATM   $(($C2+1))  P91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM    $(($C2+1))  P   UNL/HETATM    $(($C2+1))  P91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM $(($C2+1))  S   UNL/HETATM $(($C2+1))  S91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM  $(($C2+1))  S   UNL/HETATM  $(($C2+1))  S91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM   $(($C2+1))  S   UNL/HETATM   $(($C2+1))  S91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM    $(($C2+1))  S   UNL/HETATM    $(($C2+1))  S91 UNL/" $linker

line1=$(awk "/UNL/{c++} c==$(($C2+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $3 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N91 /" $linker
sed -i "${line1}s/  N${init} /  N91 /" $linker
sed -i "${line1}s/  C${init}  /  C91 /" $linker
sed -i "${line1}s/  C${init} /  C91 /" $linker
sed -i "${line1}s/  O${init}  /  O91 /" $linker
sed -i "${line1}s/  O${init} /  O91 /" $linker
sed -i "${line1}s/  P${init}  /  P91 /" $linker
sed -i "${line1}s/  P${init} /  P91 /" $linker
sed -i "${line1}s/  S${init}  /  S91 /" $linker
sed -i "${line1}s/  S${init} /  S91 /" $linker

line1=$(awk "/UNL/{c++} c==$(($C2+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $2 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N91 /" $linker
sed -i "${line1}s/  N${init} /  N91 /" $linker
sed -i "${line1}s/  C${init}  /  C91 /" $linker
sed -i "${line1}s/  C${init} /  C91 /" $linker
sed -i "${line1}s/  O${init}  /  O91 /" $linker
sed -i "${line1}s/  O${init} /  O91 /" $linker
sed -i "${line1}s/  P${init}  /  P91 /" $linker
sed -i "${line1}s/  P${init} /  P91 /" $linker
sed -i "${line1}s/  S${init}  /  S91 /" $linker
sed -i "${line1}s/  S${init} /  S91 /" $linker
fi


if [ -z "$C3" ];
then
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  N   UNL/HETATM $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  N   UNL/HETATM  $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  N   UNL/HETATM   $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  N   UNL/HETATM    $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  C   UNL/HETATM $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  C   UNL/HETATM  $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  C   UNL/HETATM   $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  C   UNL/HETATM    $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  O   UNL/HETATM $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  O   UNL/HETATM  $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  O   UNL/HETATM   $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  O   UNL/HETATM    $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  P   UNL/HETATM $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  P   UNL/HETATM  $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  P   UNL/HETATM   $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  P   UNL/HETATM    $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  S   UNL/HETATM $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  S   UNL/HETATM  $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  S   UNL/HETATM   $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  S   UNL/HETATM    $(($C5+1))  S94 UNL/" $linker

sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  N92 UNL/HETATM $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  N92 UNL/HETATM  $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  N92 UNL/HETATM   $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  N92 UNL/HETATM    $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  C92 UNL/HETATM $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  C92 UNL/HETATM  $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  C92 UNL/HETATM   $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  C92 UNL/HETATM    $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  O92 UNL/HETATM $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  O92 UNL/HETATM  $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  O92 UNL/HETATM   $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  O92 UNL/HETATM    $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  P92 UNL/HETATM $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  P92 UNL/HETATM  $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  P92 UNL/HETATM   $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  P92 UNL/HETATM    $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  S92 UNL/HETATM $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  S92 UNL/HETATM  $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  S92 UNL/HETATM   $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  S92 UNL/HETATM    $(($C5+1))  S94 UNL/" $linker

line1=$(awk "/UNL/{c++} c==$(($C5+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $3 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N94 /" $linker
sed -i "${line1}s/  N${init} /  N94 /" $linker
sed -i "${line1}s/  C${init}  /  C94 /" $linker
sed -i "${line1}s/  C${init} /  C94 /" $linker
sed -i "${line1}s/  O${init}  /  O94 /" $linker
sed -i "${line1}s/  O${init} /  O94 /" $linker
sed -i "${line1}s/  P${init}  /  P94 /" $linker
sed -i "${line1}s/  P${init} /  P94 /" $linker
sed -i "${line1}s/  S${init}  /  S94 /" $linker
sed -i "${line1}s/  S${init} /  S94 /" $linker


line1=$(awk "/UNL/{c++} c==$(($C5+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $2 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N94 /" $linker
sed -i "${line1}s/  N${init} /  N94 /" $linker
sed -i "${line1}s/  C${init}  /  C94 /" $linker
sed -i "${line1}s/  C${init} /  C94 /" $linker
sed -i "${line1}s/  O${init}  /  O94 /" $linker
sed -i "${line1}s/  O${init} /  O94 /" $linker
sed -i "${line1}s/  P${init}  /  P94 /" $linker
sed -i "${line1}s/  P${init} /  P94 /" $linker
sed -i "${line1}s/  S${init}  /  S94 /" $linker
sed -i "${line1}s/  S${init} /  S94 /" $linker
else
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  N   UNL/HETATM $(($C3+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  N   UNL/HETATM  $(($C3+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  N   UNL/HETATM   $(($C3+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  N   UNL/HETATM    $(($C3+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  C   UNL/HETATM $(($C3+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  C   UNL/HETATM  $(($C3+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  C   UNL/HETATM   $(($C3+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  C   UNL/HETATM    $(($C3+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  O   UNL/HETATM $(($C3+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  O   UNL/HETATM  $(($C3+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  O   UNL/HETATM   $(($C3+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  O   UNL/HETATM    $(($C3+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  P   UNL/HETATM $(($C3+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  P   UNL/HETATM  $(($C3+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  P   UNL/HETATM   $(($C3+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  P   UNL/HETATM    $(($C3+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  S   UNL/HETATM $(($C3+1))  S92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  S   UNL/HETATM  $(($C3+1))  S92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  S   UNL/HETATM   $(($C3+1))  S92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  S   UNL/HETATM    $(($C3+1))  S92 UNL/" $linker

line1=$(awk "/UNL/{c++} c==$(($C3+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $3 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N92 /" $linker
sed -i "${line1}s/  N${init} /  N92 /" $linker
sed -i "${line1}s/  C${init}  /  C92 /" $linker
sed -i "${line1}s/  C${init} /  C92 /" $linker
sed -i "${line1}s/  O${init}  /  O92 /" $linker
sed -i "${line1}s/  O${init} /  O92 /" $linker
sed -i "${line1}s/  P${init}  /  P92 /" $linker
sed -i "${line1}s/  P${init} /  P92 /" $linker
sed -i "${line1}s/  S${init}  /  S92 /" $linker
sed -i "${line1}s/  S${init} /  S92 /" $linker


line1=$(awk "/UNL/{c++} c==$(($C3+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $2 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N92 /" $linker
sed -i "${line1}s/  N${init} /  N92 /" $linker
sed -i "${line1}s/  C${init}  /  C92 /" $linker
sed -i "${line1}s/  C${init} /  C92 /" $linker
sed -i "${line1}s/  O${init}  /  O92 /" $linker
sed -i "${line1}s/  O${init} /  O92 /" $linker
sed -i "${line1}s/  P${init}  /  P92 /" $linker
sed -i "${line1}s/  P${init} /  P92 /" $linker
sed -i "${line1}s/  S${init}  /  S92 /" $linker
sed -i "${line1}s/  S${init} /  S92 /" $linker
fi

if [ -z "$C4" ];
then
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  N   UNL/HETATM $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  N   UNL/HETATM  $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  N   UNL/HETATM   $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  N   UNL/HETATM    $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  C   UNL/HETATM $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  C   UNL/HETATM  $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  C   UNL/HETATM   $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  C   UNL/HETATM    $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  O   UNL/HETATM $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  O   UNL/HETATM  $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  O   UNL/HETATM   $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  O   UNL/HETATM    $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  P   UNL/HETATM $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  P   UNL/HETATM  $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  P   UNL/HETATM   $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  P   UNL/HETATM    $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  S   UNL/HETATM $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  S   UNL/HETATM  $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  S   UNL/HETATM   $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  S   UNL/HETATM    $(($C5+1))  S94 UNL/" $linker

line1=$(awk "/UNL/{c++} c==$(($C5+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $3 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N94 /" $linker
sed -i "${line1}s/  N${init} /  N94 /" $linker
sed -i "${line1}s/  C${init}  /  C94 /" $linker
sed -i "${line1}s/  C${init} /  C94 /" $linker
sed -i "${line1}s/  O${init}  /  O94 /" $linker
sed -i "${line1}s/  O${init} /  O94 /" $linker
sed -i "${line1}s/  P${init}  /  P94 /" $linker
sed -i "${line1}s/  P${init} /  P94 /" $linker
sed -i "${line1}s/  S${init}  /  S94 /" $linker
sed -i "${line1}s/  S${init} /  S94 /" $linker


line1=$(awk "/UNL/{c++} c==$(($C5+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $2 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N94 /" $linker
sed -i "${line1}s/  N${init} /  N94 /" $linker
sed -i "${line1}s/  C${init}  /  C94 /" $linker
sed -i "${line1}s/  C${init} /  C94 /" $linker
sed -i "${line1}s/  O${init}  /  O94 /" $linker
sed -i "${line1}s/  O${init} /  O94 /" $linker
sed -i "${line1}s/  P${init}  /  P94 /" $linker
sed -i "${line1}s/  P${init} /  P94 /" $linker
sed -i "${line1}s/  S${init}  /  S94 /" $linker
sed -i "${line1}s/  S${init} /  S94 /" $linker
else
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM $(($C4+1))  N   UNL/HETATM $(($C4+1))  N91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM  $(($C4+1))  N   UNL/HETATM  $(($C4+1))  N91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM   $(($C4+1))  N   UNL/HETATM   $(($C4+1))  N91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM    $(($C4+1))  N   UNL/HETATM    $(($C4+1))  N91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM $(($C4+1))  C   UNL/HETATM $(($C4+1))  C91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM  $(($C4+1))  C   UNL/HETATM  $(($C4+1))  C91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM   $(($C4+1))  C   UNL/HETATM   $(($C4+1))  C91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM    $(($C4+1))  C   UNL/HETATM    $(($C4+1))  C91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM $(($C4+1))  O   UNL/HETATM $(($C4+1))  O91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM  $(($C4+1))  O   UNL/HETATM  $(($C4+1))  O91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM   $(($C4+1))  O   UNL/HETATM   $(($C4+1))  O91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM    $(($C4+1))  O   UNL/HETATM    $(($C4+1))  O91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM $(($C4+1))  P   UNL/HETATM $(($C4+1))  P91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM  $(($C4+1))  P   UNL/HETATM  $(($C4+1))  P91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM   $(($C4+1))  P   UNL/HETATM   $(($C4+1))  P91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM    $(($C4+1))  P   UNL/HETATM    $(($C4+1))  P91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM $(($C4+1))  S   UNL/HETATM $(($C4+1))  S91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM  $(($C4+1))  S   UNL/HETATM  $(($C4+1))  S91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM   $(($C4+1))  S   UNL/HETATM   $(($C4+1))  S91 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM    $(($C4+1))  S   UNL/HETATM    $(($C4+1))  S91 UNL/" $linker

line1=$(awk "/UNL/{c++} c==$(($C4+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $3 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N91 /" $linker
sed -i "${line1}s/  N${init} /  N91 /" $linker
sed -i "${line1}s/  C${init}  /  C91 /" $linker
sed -i "${line1}s/  C${init} /  C91 /" $linker
sed -i "${line1}s/  O${init}  /  O91 /" $linker
sed -i "${line1}s/  O${init} /  O91 /" $linker
sed -i "${line1}s/  P${init}  /  P91 /" $linker
sed -i "${line1}s/  P${init} /  P91 /" $linker
sed -i "${line1}s/  S${init}  /  S91 /" $linker
sed -i "${line1}s/  S${init} /  S91 /" $linker

line1=$(awk "/UNL/{c++} c==$(($C4+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $2 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N91 /" $linker
sed -i "${line1}s/  N${init} /  N91 /" $linker
sed -i "${line1}s/  C${init}  /  C91 /" $linker
sed -i "${line1}s/  C${init} /  C91 /" $linker
sed -i "${line1}s/  O${init}  /  O91 /" $linker
sed -i "${line1}s/  O${init} /  O91 /" $linker
sed -i "${line1}s/  P${init}  /  P91 /" $linker
sed -i "${line1}s/  P${init} /  P91 /" $linker
sed -i "${line1}s/  S${init}  /  S91 /" $linker
sed -i "${line1}s/  S${init} /  S91 /" $linker
fi;

if [ ! -z "$C3" ] && [ "$C2" = "$C3" ]; then

sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  N92 UNL/HETATM $(($C3+1))  N93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  N92 UNL/HETATM  $(($C3+1))  N93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  N92 UNL/HETATM   $(($C3+1))  N93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  N92 UNL/HETATM    $(($C3+1))  N93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  C92 UNL/HETATM $(($C3+1))  C93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  C92 UNL/HETATM  $(($C3+1))  C93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  C92 UNL/HETATM   $(($C3+1))  C93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  C92 UNL/HETATM    $(($C3+1))  C93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  O92 UNL/HETATM $(($C3+1))  O93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  O92 UNL/HETATM  $(($C3+1))  O93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  O92 UNL/HETATM   $(($C3+1))  O93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  O92 UNL/HETATM    $(($C3+1))  O93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  P92 UNL/HETATM $(($C3+1))  P93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  P92 UNL/HETATM  $(($C3+1))  P93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  P92 UNL/HETATM   $(($C3+1))  P93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  P92 UNL/HETATM    $(($C3+1))  P93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  S92 UNL/HETATM $(($C3+1))  S93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  S92 UNL/HETATM  $(($C3+1))  S93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  S92 UNL/HETATM   $(($C3+1))  S93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  S92 UNL/HETATM    $(($C3+1))  S93 UNL/" $linker

line1=$(awk "/UNL/{c++} c==$(($C3+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $3 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N93 /" $linker
sed -i "${line1}s/  N${init} /  N93 /" $linker
sed -i "${line1}s/  C${init}  /  C93 /" $linker
sed -i "${line1}s/  C${init} /  C93 /" $linker
sed -i "${line1}s/  O${init}  /  O93 /" $linker
sed -i "${line1}s/  O${init} /  O93 /" $linker
sed -i "${line1}s/  P${init}  /  P93 /" $linker
sed -i "${line1}s/  P${init} /  P93 /" $linker
sed -i "${line1}s/  S${init}  /  S93 /" $linker
sed -i "${line1}s/  S${init} /  S93 /" $linker


line1=$(awk "/UNL/{c++} c==$(($C3+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $2 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N93 /" $linker
sed -i "${line1}s/  N${init} /  N93 /" $linker
sed -i "${line1}s/  C${init}  /  C93 /" $linker
sed -i "${line1}s/  C${init} /  C93 /" $linker
sed -i "${line1}s/  O${init}  /  O93 /" $linker
sed -i "${line1}s/  O${init} /  O93 /" $linker
sed -i "${line1}s/  P${init}  /  P93 /" $linker
sed -i "${line1}s/  P${init} /  P93 /" $linker
sed -i "${line1}s/  S${init}  /  S93 /" $linker
sed -i "${line1}s/  S${init} /  S93 /" $linker


sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM $(($C4+1))  N91 UNL/HETATM $(($C4+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM  $(($C4+1))  N91 UNL/HETATM  $(($C4+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM   $(($C4+1))  N91 UNL/HETATM   $(($C4+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM    $(($C4+1))  N91 UNL/HETATM    $(($C4+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM $(($C4+1))  C91 UNL/HETATM $(($C4+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM  $(($C4+1))  C91 UNL/HETATM  $(($C4+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM   $(($C4+1))  C91 UNL/HETATM   $(($C4+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM    $(($C4+1))  C91 UNL/HETATM    $(($C4+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM $(($C4+1))  O91 UNL/HETATM $(($C4+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM  $(($C4+1))  O91 UNL/HETATM  $(($C4+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM   $(($C4+1))  O91 UNL/HETATM   $(($C4+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM    $(($C4+1))  O91 UNL/HETATM    $(($C4+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM $(($C4+1))  P91 UNL/HETATM $(($C4+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM  $(($C4+1))  P91 UNL/HETATM  $(($C4+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM   $(($C4+1))  P91 UNL/HETATM   $(($C4+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM    $(($C4+1))  P91 UNL/HETATM    $(($C4+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM $(($C4+1))  S91 UNL/HETATM $(($C4+1))  S92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM  $(($C4+1))  S91 UNL/HETATM  $(($C4+1))  S92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM   $(($C4+1))  S91 UNL/HETATM   $(($C4+1))  S92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C4+1))" -e "s/HETATM    $(($C4+1))  S91 UNL/HETATM    $(($C4+1))  S92 UNL/" $linker

line1=$(awk "/UNL/{c++} c==$(($C4+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $3 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N92 /" $linker
sed -i "${line1}s/  N${init} /  N92 /" $linker
sed -i "${line1}s/  C${init}  /  C92 /" $linker
sed -i "${line1}s/  C${init} /  C92 /" $linker
sed -i "${line1}s/  O${init}  /  O92 /" $linker
sed -i "${line1}s/  O${init} /  O92 /" $linker
sed -i "${line1}s/  P${init}  /  P92 /" $linker
sed -i "${line1}s/  P${init} /  P92 /" $linker
sed -i "${line1}s/  S${init}  /  S92 /" $linker
sed -i "${line1}s/  S${init} /  S92 /" $linker


line1=$(awk "/UNL/{c++} c==$(($C4+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $2 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N92 /" $linker
sed -i "${line1}s/  N${init} /  N92 /" $linker
sed -i "${line1}s/  C${init}  /  C92 /" $linker
sed -i "${line1}s/  C${init} /  C92 /" $linker
sed -i "${line1}s/  O${init}  /  O92 /" $linker
sed -i "${line1}s/  O${init} /  O92 /" $linker
sed -i "${line1}s/  P${init}  /  P92 /" $linker
sed -i "${line1}s/  P${init} /  P92 /" $linker
sed -i "${line1}s/  S${init}  /  S92 /" $linker
sed -i "${line1}s/  S${init} /  S92 /" $linker
else
:
fi;

if [ ! -z "$C4" ] && [ "$C3" = "$C4" ]; then

sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  N92 UNL/HETATM $(($C3+1))  N93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  N92 UNL/HETATM  $(($C3+1))  N93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  N92 UNL/HETATM   $(($C3+1))  N93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  N92 UNL/HETATM    $(($C3+1))  N93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  C92 UNL/HETATM $(($C3+1))  C93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  C92 UNL/HETATM  $(($C3+1))  C93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  C92 UNL/HETATM   $(($C3+1))  C93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  C92 UNL/HETATM    $(($C3+1))  C93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  O92 UNL/HETATM $(($C3+1))  O93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  O92 UNL/HETATM  $(($C3+1))  O93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  O92 UNL/HETATM   $(($C3+1))  O93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  O92 UNL/HETATM    $(($C3+1))  O93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  P92 UNL/HETATM $(($C3+1))  P93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  P92 UNL/HETATM  $(($C3+1))  P93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  P92 UNL/HETATM   $(($C3+1))  P93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  P92 UNL/HETATM    $(($C3+1))  P93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM $(($C3+1))  S92 UNL/HETATM $(($C3+1))  S93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM  $(($C3+1))  S92 UNL/HETATM  $(($C3+1))  S93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM   $(($C3+1))  S92 UNL/HETATM   $(($C3+1))  S93 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C3+1))" -e "s/HETATM    $(($C3+1))  S92 UNL/HETATM    $(($C3+1))  S93 UNL/" $linker

line1=$(awk "/UNL/{c++} c==$(($C3+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $3 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N93 /" $linker
sed -i "${line1}s/  N${init} /  N93 /" $linker
sed -i "${line1}s/  C${init}  /  C93 /" $linker
sed -i "${line1}s/  C${init} /  C93 /" $linker
sed -i "${line1}s/  O${init}  /  O93 /" $linker
sed -i "${line1}s/  O${init} /  O93 /" $linker
sed -i "${line1}s/  P${init}  /  P93 /" $linker
sed -i "${line1}s/  P${init} /  P93 /" $linker
sed -i "${line1}s/  S${init}  /  S93 /" $linker
sed -i "${line1}s/  S${init} /  S93 /" $linker



line1=$(awk "/UNL/{c++} c==$(($C3+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $2 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N93 /" $linker
sed -i "${line1}s/  N${init} /  N93 /" $linker
sed -i "${line1}s/  C${init}  /  C93 /" $linker
sed -i "${line1}s/  C${init} /  C93 /" $linker
sed -i "${line1}s/  O${init}  /  O93 /" $linker
sed -i "${line1}s/  O${init} /  O93 /" $linker
sed -i "${line1}s/  P${init}  /  P93 /" $linker
sed -i "${line1}s/  P${init} /  P93 /" $linker
sed -i "${line1}s/  S${init}  /  S93 /" $linker
sed -i "${line1}s/  S${init} /  S93 /" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM $(($C2+1))  N91 UNL/HETATM $(($C2+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM  $(($C2+1))  N91 UNL/HETATM  $(($C2+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM   $(($C2+1))  N91 UNL/HETATM   $(($C2+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM    $(($C2+1))  N91 UNL/HETATM    $(($C2+1))  N92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM $(($C2+1))  C91 UNL/HETATM $(($C2+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM  $(($C2+1))  C91 UNL/HETATM  $(($C2+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM   $(($C2+1))  C91 UNL/HETATM   $(($C2+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM    $(($C2+1))  C91 UNL/HETATM    $(($C2+1))  C92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM $(($C2+1))  O91 UNL/HETATM $(($C2+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM  $(($C2+1))  O91 UNL/HETATM  $(($C2+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM   $(($C2+1))  O91 UNL/HETATM   $(($C2+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM    $(($C2+1))  O91 UNL/HETATM    $(($C2+1))  O92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM $(($C2+1))  P91 UNL/HETATM $(($C2+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM  $(($C2+1))  P91 UNL/HETATM  $(($C2+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM   $(($C2+1))  P91 UNL/HETATM   $(($C2+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM    $(($C2+1))  P91 UNL/HETATM    $(($C2+1))  P92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM $(($C2+1))  S91 UNL/HETATM $(($C2+1))  S92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM  $(($C2+1))  S91 UNL/HETATM  $(($C2+1))  S92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM   $(($C2+1))  S91 UNL/HETATM   $(($C2+1))  S92 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C2+1))" -e "s/HETATM    $(($C2+1))  S91 UNL/HETATM    $(($C2+1))  S92 UNL/" $linker

line1=$(awk "/UNL/{c++} c==$(($C2+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $3 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N92 /" $linker
sed -i "${line1}s/  N${init} /  N92 /" $linker
sed -i "${line1}s/  C${init}  /  C92 /" $linker
sed -i "${line1}s/  C${init} /  C92 /" $linker
sed -i "${line1}s/  O${init}  /  O92 /" $linker
sed -i "${line1}s/  O${init} /  O92 /" $linker
sed -i "${line1}s/  P${init}  /  P92 /" $linker
sed -i "${line1}s/  P${init} /  P92 /" $linker
sed -i "${line1}s/  S${init}  /  S92 /" $linker
sed -i "${line1}s/  S${init} /  S92 /" $linker

line1=$(awk "/UNL/{c++} c==$(($C2+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $2 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N92 /" $linker
sed -i "${line1}s/  N${init} /  N92 /" $linker
sed -i "${line1}s/  C${init}  /  C92 /" $linker
sed -i "${line1}s/  C${init} /  C92 /" $linker
sed -i "${line1}s/  O${init}  /  O92 /" $linker
sed -i "${line1}s/  O${init} /  O92 /" $linker
sed -i "${line1}s/  P${init}  /  P92 /" $linker
sed -i "${line1}s/  P${init} /  P92 /" $linker
sed -i "${line1}s/  S${init}  /  S92 /" $linker
sed -i "${line1}s/  S${init} /  S92 /" $linker
else
:
fi;
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  N   UNL/HETATM $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  N   UNL/HETATM  $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  N   UNL/HETATM   $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  N   UNL/HETATM    $(($C5+1))  N94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  C   UNL/HETATM $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  C   UNL/HETATM  $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  C   UNL/HETATM   $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  C   UNL/HETATM    $(($C5+1))  C94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  O   UNL/HETATM $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  O   UNL/HETATM  $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  O   UNL/HETATM   $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  O   UNL/HETATM    $(($C5+1))  O94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  P   UNL/HETATM $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  P   UNL/HETATM  $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  P   UNL/HETATM   $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  P   UNL/HETATM    $(($C5+1))  P94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM $(($C5+1))  S   UNL/HETATM $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM  $(($C5+1))  S   UNL/HETATM  $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM   $(($C5+1))  S   UNL/HETATM   $(($C5+1))  S94 UNL/" $linker
sed -i -z -e "s/UNL/UNL/$(($C5+1))" -e "s/HETATM    $(($C5+1))  S   UNL/HETATM    $(($C5+1))  S94 UNL/" $linker

line1=$(awk "/UNL/{c++} c==$(($C5+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $3 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N94 /" $linker
sed -i "${line1}s/  N${init} /  N94 /" $linker
sed -i "${line1}s/  C${init}  /  C94 /" $linker
sed -i "${line1}s/  C${init} /  C94 /" $linker
sed -i "${line1}s/  O${init}  /  O94 /" $linker
sed -i "${line1}s/  O${init} /  O94 /" $linker
sed -i "${line1}s/  P${init}  /  P94 /" $linker
sed -i "${line1}s/  P${init} /  P94 /" $linker
sed -i "${line1}s/  S${init}  /  S94 /" $linker
sed -i "${line1}s/  S${init} /  S94 /" $linker


line1=$(awk "/UNL/{c++} c==$(($C5+1)){print NR;exit}" $linker)
C9=$(awk "NR==$line1 { print $2 }" $linker)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N94 /" $linker
sed -i "${line1}s/  N${init} /  N94 /" $linker
sed -i "${line1}s/  C${init}  /  C94 /" $linker
sed -i "${line1}s/  C${init} /  C94 /" $linker
sed -i "${line1}s/  O${init}  /  O94 /" $linker
sed -i "${line1}s/  O${init} /  O94 /" $linker
sed -i "${line1}s/  P${init}  /  P94 /" $linker
sed -i "${line1}s/  P${init} /  P94 /" $linker
sed -i "${line1}s/  S${init}  /  S94 /" $linker
sed -i "${line1}s/  S${init} /  S94 /" $linker
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM $(($C6+1))  N   UNL/HETATM $(($C6+1))  N95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM  $(($C6+1))  N   UNL/HETATM  $(($C6+1))  N95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM   $(($C6+1))  N   UNL/HETATM   $(($C6+1))  N95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM    $(($C6+1))  N   UNL/HETATM    $(($C6+1))  N95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM $(($C6+1))  C   UNL/HETATM $(($C6+1))  C95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM  $(($C6+1))  C   UNL/HETATM  $(($C6+1))  C95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM   $(($C6+1))  C   UNL/HETATM   $(($C6+1))  C95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM    $(($C6+1))  C   UNL/HETATM    $(($C6+1))  C95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM $(($C6+1))  O   UNL/HETATM $(($C6+1))  O95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM  $(($C6+1))  O   UNL/HETATM  $(($C6+1))  O95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM   $(($C6+1))  O   UNL/HETATM   $(($C6+1))  O95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM    $(($C6+1))  O   UNL/HETATM    $(($C6+1))  O95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM $(($C6+1))  P   UNL/HETATM $(($C6+1))  P95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM  $(($C6+1))  P   UNL/HETATM  $(($C6+1))  P95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM   $(($C6+1))  P   UNL/HETATM   $(($C6+1))  P95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM    $(($C6+1))  P   UNL/HETATM    $(($C6+1))  P95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM $(($C6+1))  S   UNL/HETATM $(($C6+1))  S95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM  $(($C6+1))  S   UNL/HETATM  $(($C6+1))  S95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM   $(($C6+1))  S   UNL/HETATM   $(($C6+1))  S95 UNL/" $e3
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM    $(($C6+1))  S   UNL/HETATM    $(($C6+1))  S95 UNL/" $e3

line1=$(awk "/UNL/{c++} c==$(($C6+1)){print NR;exit}" $e3)
C9=$(awk "NR==$line1 { print $3 }" $e3)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N95 /" $e3
sed -i "${line1}s/  N${init} /  N95 /" $e3
sed -i "${line1}s/  C${init}  /  C95 /" $e3
sed -i "${line1}s/  C${init} /  C95 /" $e3
sed -i "${line1}s/  O${init}  /  O95 /" $e3
sed -i "${line1}s/  O${init} /  O95 /" $e3
sed -i "${line1}s/  P${init}  /  P95 /" $e3
sed -i "${line1}s/  P${init} /  P95 /" $e3
sed -i "${line1}s/  S${init}  /  S95 /" $e3
sed -i "${line1}s/  S${init} /  S95 /" $e3

line1=$(awk "/UNL/{c++} c==$(($C6+1)){print NR;exit}" $e3)
C9=$(awk "NR==$line1 { print $2 }" $e3)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N95 /" $e3
sed -i "${line1}s/  N${init} /  N95 /" $e3
sed -i "${line1}s/  C${init}  /  C95 /" $e3
sed -i "${line1}s/  C${init} /  C95 /" $e3
sed -i "${line1}s/  O${init}  /  O95 /" $e3
sed -i "${line1}s/  O${init} /  O95 /" $e3
sed -i "${line1}s/  P${init}  /  P95 /" $e3
sed -i "${line1}s/  P${init} /  P95 /" $e3
sed -i "${line1}s/  S${init}  /  S95 /" $e3
sed -i "${line1}s/  S${init} /  S95 /" $e3

done;
