#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=LEAP

source ~/.bashrc
conda activate amber


p=$1
e=$2


cp ${p%.*}_r_prep_pp_p1fixed_filtered.mol2 ${p%.*}_r_prep_pp_fixed.mol2
cp ${e%.*}_l_prep_pp_p1fixed_filtered.mol2 ${e%.*}_l_prep_pp_fixed.mol2

poi="${p%.*}_r_prep_pp_fixed.mol2"
e3="${e%.*}_l_prep_pp_fixed.mol2"

rm xx* *model* diagnosis* *temp*
csplit $e3 '/@<TRIPOS>MOLECULE/' '{*}'
rm xx00 

#poi_7pi4_docked_rec_r_prep_pp_p1fixed_filtered.mol2

sed -i -e '/UNL/ s/ATOM  /HETATM/'
obabel -imol2 $poi -opdb -O ${poi%.*}.pdb
sed -i -e '/UNL/ s/ATOM  /HETATM/' ${poi%.*}.pdb
sed -i '/HETATM/d' ${poi%.*}.pdb
sed -i '/CONECT/d' ${poi%.*}.pdb
sed -i '/MASTER/d' ${poi%.*}.pdb
sed -i '/END/d' ${poi%.*}.pdb

v1=$(echo "${poi%.*}.pdb" | cut -d '_' -f 2)
v2=$(echo "${poi%.*}.pdb" | cut -d '_' -f 3)
v30=$(echo "$e3" | cut -d '.' -f1 | rev | cut -d '_' -f7-8 | rev)
v31=$(echo "$poi" | cut -d '.' -f1 | rev | cut -d '_' -f5-6 | rev)
v3=${v30}_${v31}

grep "^HETATM" $p > ${poi%.*}_temp_lig.pdb
grep "^HETATM" $e > ${v3}_temp_lig2.pdb
grep "^CONECT" $e >> ${v3}_temp_lig2.pdb
obabel -ipdb ${v3}_temp_lig2.pdb -osd -O ${v3}_temp_lig2.sdf

for x in xx*; do
obabel -imol2 $x -opdb -O ${x%.*}.pdb
sed -i -e '/UNL/ s/ATOM  /HETATM/' ${x%.*}.pdb
sed -i '/ UNL /{s/ A /   /}' ${x%.*}.pdb

cp ${x%.*}.pdb Te3G_${x%.*}.pdb
pdb_base="Te3G_${x%.*}"
up_pdb=${pdb_base^^}
ligand="unl"
$ICMHOME/icm -vlscluster $ICMHOME/_dockBatch Te3G_${x%.*}.pdb ligmol=$ligand;
rm -r ${pdb_base}_maps/;
mkdir ${pdb_base}_maps/;
mv -t ${pdb_base}_maps/ *${up_pdb}_*map *${up_pdb}_*htm *${up_pdb}.dtb *${up_pdb}_*ob;
echo "initiating docking for .." $ligand;

###VARIABLES
project=${pdb_base}_maps/D_${up_pdb}
input_sdf=${v3}_temp_lig2.sdf
thoroughness=8.0
output=${v3}_temp_lig2_docked
size_single=2                              #compounds in each job (first job will have $size_single-1)

###ICM RUN
$ICMHOME/icm -vlscluster $ICMHOME/_dockScan $project -a input=$input_sdf maxFileSizeMb=5000 effort=$thoroughness name=$output output=$output >> $output.log
rm -r ${pdb_base}_maps/
obabel -isd ${v3}_temp_lig2_docked.sdf -opdb -O ${v3}_temp_lig2_docked.pdb
	
sed -i '/CONECT/d' ${v3}_temp_lig2_docked.pdb
sed -i '/COMPND/d' ${v3}_temp_lig2_docked.pdb
sed -i '/AUTHOR/d' ${v3}_temp_lig2_docked.pdb
sed -i '/MASTER/d' ${v3}_temp_lig2_docked.pdb
sed -i '/END/d' ${v3}_temp_lig2_docked.pdb

up_ligand="UNL"
sed -i "/HETATM/{/${up_ligand}/d}" ${x%.*}.pdb
sed -i '/CONECT/d' ${x%.*}.pdb 
sed -i '/MASTER/d' ${x%.*}.pdb
sed -i '/END/d' ${x%.*}.pdb
sed -i '/COMPND/d' ${x%.*}.pdb
sed -i '/AUTHOR/d' ${x%.*}.pdb

v1l=$(echo "$v30" | tr '[:lower:]' '[:upper:]')

con=$(ls ../conect_cp/${v1l}_cnct.idx)

## amending pdb for linker ligation ## LINKER FIRST ##
line=$(head -n 1 $con)
C1=$(cut -d',' -f1 <<<$line)
C2=$(cut -d',' -f2 <<<$line)
C3=$(cut -d',' -f3 <<<$line)
C4=$(cut -d',' -f4 <<<$line)
C5=$(cut -d',' -f5 <<<$line)
C6=$(cut -d',' -f6 <<<$line)

sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM $(($C6+1))  N   UNL/HETATM $(($C6+1))  N95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM  $(($C6+1))  N   UNL/HETATM  $(($C6+1))  N95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM   $(($C6+1))  N   UNL/HETATM   $(($C6+1))  N95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM    $(($C6+1))  N   UNL/HETATM    $(($C6+1))  N95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM $(($C6+1))  C   UNL/HETATM $(($C6+1))  C95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM  $(($C6+1))  C   UNL/HETATM  $(($C6+1))  C95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM   $(($C6+1))  C   UNL/HETATM   $(($C6+1))  C95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM    $(($C6+1))  C   UNL/HETATM    $(($C6+1))  C95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM $(($C6+1))  O   UNL/HETATM $(($C6+1))  O95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM  $(($C6+1))  O   UNL/HETATM  $(($C6+1))  O95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM   $(($C6+1))  O   UNL/HETATM   $(($C6+1))  O95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM    $(($C6+1))  O   UNL/HETATM    $(($C6+1))  O95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM $(($C6+1))  P   UNL/HETATM $(($C6+1))  P95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM  $(($C6+1))  P   UNL/HETATM  $(($C6+1))  P95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM   $(($C6+1))  P   UNL/HETATM   $(($C6+1))  P95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM    $(($C6+1))  P   UNL/HETATM    $(($C6+1))  P95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM $(($C6+1))  S   UNL/HETATM $(($C6+1))  S95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM  $(($C6+1))  S   UNL/HETATM  $(($C6+1))  S95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM   $(($C6+1))  S   UNL/HETATM   $(($C6+1))  S95 UNL/" ${v3}_temp_lig2_docked.pdb
sed -i -z -e "s/UNL/UNL/$(($C6+1))" -e "s/HETATM    $(($C6+1))  S   UNL/HETATM    $(($C6+1))  S95 UNL/" ${v3}_temp_lig2_docked.pdb

line1=$(awk "/UNL/{c++} c==$(($C6+1)){print NR;exit}" ${v3}_temp_lig2_docked.pdb)

C9=$(awk "NR==$line1 { print $3 }" ${v3}_temp_lig2_docked.pdb)
C9b=$(echo $C9 | cut -d ' ' -f3)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  N${init} /  N95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  C${init}  /  C95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  C${init} /  C95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  O${init}  /  O95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  O${init} /  O95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  P${init}  /  P95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  P${init} /  P95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  S${init}  /  S95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  S${init} /  S95 /" ${v3}_temp_lig2_docked.pdb

line1=$(awk "/UNL/{c++} c==$(($C6+1)){print NR;exit}" ${v3}_temp_lig2_docked.pdb)
C9=$(awk "NR==$line1 { print $2 }" ${v3}_temp_lig2_docked.pdb)
C9b=$(echo $C9 | cut -d ' ' -f2)
init=${C9b:1:3}
sed -i "${line1}s/  N${init}  /  N95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  N${init} /  N95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  C${init}  /  C95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  C${init} /  C95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  O${init}  /  O95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  O${init} /  O95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  P${init}  /  P95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  P${init} /  P95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  S${init}  /  S95 /" ${v3}_temp_lig2_docked.pdb
sed -i "${line1}s/  S${init} /  S95 /" ${v3}_temp_lig2_docked.pdb


cat ${poi%.*}.pdb ${poi%.*}_temp_lig.pdb ${x%.*}.pdb ${v3}_temp_lig2_docked.pdb >> complex_${v3}_pp_model_${x%.*}.pdb
pdb4amber -i complex_${v3}_pp_model_${x%.*}.pdb -o amb_complex_${v3}_pp_model_${x%.*}.pdb
mv amb_complex_${v3}_pp_model_${x%.*}.pdb complex_${v3}_pp_model_${x%.*}.pdb
cp complex_${v3}_pp_model_${x%.*}.pdb ../top_20_pooled

rm -r  ${v3}_temp_lig2_docked.sdf *maps* ${v3}_temp_lig2_docked.pdb ${x%.*}.pdb *log *XX* *amb*

done;

rm ${poi%.*}_temp_lig.pdb *xx* *temp* 



cd ../top_20_pooled
for c in complex_${v3}_pp_model*pdb; do sbatch maestrov2.sh $c;done

