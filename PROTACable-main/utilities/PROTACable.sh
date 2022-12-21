#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=protacs

source ~/.bashrc
conda activate py39


###### INSTRUCTIONS ######
# 1) Put all .mol fragment in one directory. label them short two_syllabyl names with meaningful chemotype description in the 2nd part.
# 2) Put the poi in that directory. Label it short name like poi_target.pdb. It must have a ligand bound to the site of "nukking".
# 3) Copy the utilities to this directory. Copy PROTACable.sh out from the utilities to this directory!
# 4) Command example: bash PROTACable.sh original_ligand_with_untampered_chemotypes.mol poi_target.pdb original_fragment_exit_vector_idx->int 2nd_fragment_amide.mol its_exit_vector->int 3rd_fragment_carboxy its_exit_vector->int



ligandPoi=$1
poi=$2
unique_idx=$3
fragt=$4
fragt_idx=$5
fragth=$6
fragth_idx=$7

genesis=$PWD
mkdir ${1%.*}_${2%.*}_${3%.*}


mkdir previous_results
mkdir previous_results/chosen_e3
mv chosen_e3/* previous_results/chosen_e3

cp ./{$1,$2} ${1%.*}_${2%.*}_${3%.*}
cd ${1%.*}_${2%.*}_${3%.*}
first_iteration="${genesis}/${1%.*}_${2%.*}_${3%.*}"

echo "genesis directory is ${genesis}"
echo "1st_iteration directory is ${first_iteration}"

rm -r -f *rms *mm *log *mol2 *sdf *maps* sed* *cp* linker_library/temp/
cp -r ../utilities/* ./
cp -r ../utilities/ ./
cp $poi ${poi%.*}_cp.pdb

pdb="${poi%.*}_cp.pdb"

## 1) dock ligand to POI

pdb_base=${pdb%.*};
grep "^HETATM" $pdb > ${pdb_base}_ligand.pdb;
ligand=$(awk 'NR==1 {print $4}' ${pdb_base}_ligand.pdb | awk '{print $0}');
sed -i "/HETATM/s/${ligand}/UNL/" $pdb
rm ${pdb_base}_ligand.pdb;
lett="A B C D E F G H"; for L in $lett; do echo $L; sed -i "/ UNL /{s/ UNL ${L} / UNL   /}" $pdb;done

$ICMHOME/icm -vlscluster $ICMHOME/_dockBatch $pdb ligmol=unl
mkdir ${pdb_base}_maps/;
mv -t ${pdb_base}_maps/ *map *htm *dtb *ob;

echo "initiating docking for .." $ligandPoi;
up_pdb=${pdb_base^^};

sbatch --wait icm_poi.sh ${pdb_base}_maps/D_${up_pdb} $ligandPoi ${ligandPoi%.*}_docked; 

for i in ${ligandPoi%.*}_docked*sdf; do mv $i ${ligandPoi%.*}_docked.sdf;done

## 2) add ligand to receptor

# remove existing ligand
grep "^HETATM" $pdb > ${pdb_base}_ligand.pdb;
#ligand=$(awk 'NR==1 {print $4}' ${pdb_base}_ligand.pdb | awk '{print tolower($0)}');
rm ${pdb_base}_ligand.pdb;
ligand="UNL"
up_ligand=${ligand^^};
sed -i "/HETATM/{/${up_ligand}/d}" $pdb;
sed -i "/HETATM/{/${ligand}/d}" $pdb;


# add new ligand
cp $pdb ${pdb_base}_rec.pdb;
sed -i '/CONECT/d' ${pdb_base}_rec.pdb
sed -i '/END/d' ${pdb_base}_rec.pdb
obabel -isd ${ligandPoi%.*}_docked.sdf -opdb -O ${ligandPoi%.*}_docked.pdb
cat ${ligandPoi%.*}_docked.pdb >> ${pdb_base}_rec.pdb;
sed -i '/AUTHOR/d' ${pdb_base}_rec.pdb
sed -i '/COMPND/d' ${pdb_base}_rec.pdb

## 3) dock linker on top of the ligand

obabel -ipdb ${ligandPoi%.*}_docked.pdb -omol2 -O ${ligandPoi%.*}_docked.mol2
        
mol2_lig=${ligandPoi}_docked.mol2
rm -r *_maps
$ICMHOME/icm -vlscluster $ICMHOME/_dockBatch ${pdb_base}_rec.pdb ligmol=${mol2_lig};
mkdir ${pdb_base}_maps/;
mv -t ${pdb_base}_maps/ *map *htm *dtb *ob;

linker_count=$(ls linker_library/* | wc -l )

#for idx in $(seq 1 $linker_count); do 
echo "DOCKING DUMMY LINKERS"

for mol in dummy_linkers/*mol; do
	echo "initiating docking for .. DUMMY linker ${mol}"
	sbatch icm_linker.sh ${pdb_base}_maps/D_${up_pdb}_REC $mol ${mol%.*}_linker_docked
done;

i=1
for mol in linker_library/*mol; do
	echo "initiating docking for .. linker ${i}";
	sbatch icm_linker.sh ${pdb_base}_maps/D_${up_pdb}_REC $mol ${mol%.*}_linker_docked
	sleep 1s
let "i+=1"
	
done;

while [[ $(squeue --name icm | wc -l) -gt 1 ]]; do echo "wait"; sleep 1s; done;

sbatch wait.sh

echo "renaming linkers"
for docked in linker_library/*sdf; do mv "$docked" "${docked%_linker_docked*.sdf}_linker_docked.sdf"; done

for docked in dummy_linkers/*sdf; do mv "$docked" "${docked%_linker_docked*.sdf}_linker_docked.sdf"; done


rm linker_library/*log
rm dummy_linkers/*log

mkdir dummy_linkers/temp
mv dummy_linkers/*sdf dummy_linkers/temp

mkdir linker_library/temp
mv linker_library/*sdf linker_library/temp


## 4) modify conect records for new POI

cp -r conect/ conect_cp/

for file in conect_cp/*idx; do
	line=$(head -n 1 $file)
	C1=$(cut -d',' -f1 <<<$line)
	C2=$(cut -d',' -f2 <<<$line)
	C3=$(cut -d',' -f3 <<<$line)
	C4=$(cut -d',' -f4 <<<$line)
	C5=$(cut -d',' -f5 <<<$line)
	C6=$(cut -d',' -f6 <<<$line)
	original="$C1,$C2,$C3,$C4,$C5,$C6"
	C1=${unique_idx}
	if [ -z "$C2" ]; then
		C5=${unique_idx}	
	else
		:
	fi
	new="$C1,$C2,$C3,$C4,$C5,$C6"
	sed -i "s/${original}/${new}/" ${file}
done;

cp -r dummy_conect/ dummy_conect_cp/

for file in dummy_conect_cp/*idx; do
        line=$(head -n 1 $file)
        C1=$(cut -d',' -f1 <<<$line)
        C2=$(cut -d',' -f2 <<<$line)
        C3=$(cut -d',' -f3 <<<$line)
        C4=$(cut -d',' -f4 <<<$line)
        C5=$(cut -d',' -f5 <<<$line)
        C6=$(cut -d',' -f6 <<<$line)
        original="$C1,$C2,$C3,$C4,$C5,$C6"
        C1=${unique_idx}
        if [ -z "$C2" ]; then
                C5=${unique_idx}
        else
                :
        fi
        new="$C1,$C2,$C3,$C4,$C5,$C6"
        sed -i "s/${original}/${new}/" ${file}
done;

## 5) select optimal distance linker 

python SelectConformer.py

## 6) stage1 minimization -- adding linkers and conect paramterization

# convert sdf linkers to pdb

for linker in linker_library/temp/*short.sdf; do obabel -isd $linker -opdb -O ${linker%.*}.pdb;done
for linker in dummy_linkers/temp/*short.sdf; do obabel -isd $linker -opdb -O ${linker%.*}.pdb;done

poi=$(ls *_cp_rec.pdb)

for linker in linker_library/temp/*short.pdb; do 


link=${linker##*/}

v1=$(echo "$link" | cut -d '_' -f 1)
v2=$(echo "$link" | cut -d '_' -f 2)
v3=${v1}_${v2}

con=$(ls conect_cp/${v3}_*idx)


## amending pdb for linker ligation ## LINKER FIRST ##
line=$(head -n 1 $con)
C1=$(cut -d',' -f1 <<<$line)
C2=$(cut -d',' -f2 <<<$line)
C3=$(cut -d',' -f3 <<<$line)
C4=$(cut -d',' -f4 <<<$line)
C5=$(cut -d',' -f5 <<<$line)
C6=$(cut -d',' -f6 <<<$line)

echo "${poi},${linker},${con},${C1},${C2},${C3},${C4},${C5},${C6}"
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


sbatch stage1.sh $poi ${PWD}/${linker}
done;
while [[ $(squeue --name minimize1 | wc -l) -gt 1 ]]; do echo "waiting for 1st stage minimization to finish"; sleep 10s; done;


for linker in dummy_linkers/temp/*short.pdb; do


link=${linker##*/}

v1=$(echo "$link" | cut -d '_' -f 1)
v2=$(echo "$link" | cut -d '_' -f 2)
v3=${v1}_${v2}

con=$(ls dummy_conect_cp/${v3}_*idx)


## amending pdb for linker ligation ## LINKER FIRST ##
line=$(head -n 1 $con)
C1=$(cut -d',' -f1 <<<$line)
C2=$(cut -d',' -f2 <<<$line)
C3=$(cut -d',' -f3 <<<$line)
C4=$(cut -d',' -f4 <<<$line)
C5=$(cut -d',' -f5 <<<$line)
C6=$(cut -d',' -f6 <<<$line)

echo "${poi},${linker},${con},${C1},${C2},${C3},${C4},${C5},${C6}"
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


sbatch stage1_dummy.sh $poi ${PWD}/${linker}
done;
while [[ $(squeue --name minimize1 | wc -l) -gt 1 ]]; do echo "waiting for 1st stage minimization for linkers to finish"; sleep 10s; done;

rm -r *mol2 *maps*
echo "ca y est"

poi=$(ls *_cp_rec.pdb)

rm -r selected1_sorted_filtered_all_rmsd.rms selected1_sorted_all_rmsd.rms sorted_all_rmsd.rms all_mmgbsa_scores_pt1.mm sorted_all_mmgbsa_scores_pt1.mm selected2_mmgbsa_pt1.mm *_mmgbsa_rescored.pdb *-out-out.maegz *.mae *out.csv *_pair
for r in *rms; do f1=$(head -n 1 $r | cut -d ' ' -f3); if [ ! -z "$f1" ]; then echo "${r} ${f1}" >>all_rmsd.rms; else : ; fi; done; 
for r in *rmsdum; do f1=$(head -n 1 $r | cut -d ' ' -f3); if [ ! -z "$f1" ]; then echo "${r} ${f1}" >>all_rmsd.rmsdum; else : ; fi; done;


sed -i '/inf/d' all_rmsd.rms; 

sort -k3 -n all_rmsd.rms >> sorted_all_rmsd.rms;
sort -k3 -n all_rmsd.rmsdum >> sorted_all_rmsd.rmsdum;

head -n 250 sorted_all_rmsd.rms| shuf | head -n 100 >> selected1_sorted_all_rmsd.rms; 

awk '{print $1,$2}' selected1_sorted_all_rmsd.rms >> selected1_sorted_filtered_all_rmsd.rms;
awk '{print $1,$2}' sorted_all_rmsd.rmsdum >> selected1_sorted_filtered_all_rmsd.rmsdum;
#sed -i -e 's/\s\+//g' selected1_sorted_all_rmsd.rms;

#while read line; do id=$(echo "$line" | cut -d ' ' -f1 | cut -d '.' -f1 | rev | cut -d '_' -f5-6 | rev); $SCHRODINGER/utilities/structconvert -ipdb out_nhyd99_mae_amb_*${id}*pdb -omae out_nhyd99_mae_amb_${poi%.*}_${id}.mae; $SCHRODINGER/run pv_convert.py out_nhyd99_mae_amb_${poi%.*}_${id}.mae -mode split_pv -lig_last_mol; $SCHRODINGER/prime_mmgbsa out_nhyd99_mae_amb_${poi%.*}_${id%.*}-out_pv.mae -HOST slurm-compute -WAIT -NJOBS 1 -report_prime_log yes -D -job_type SITE_OPT; $SCHRODINGER/utilities/structconvert -imae out_nhyd99_mae_amb_${poi%.*}_${id%.*}-out-out.maegz -opdb poi_${poi%.*}_${id%.*}_mmgbsa_rescored.pdb; done <selected1_sorted_filtered_all_rmsd.rms;

rm *_pair/*out
#for r in *-out-out.csv; do f1=$(awk -F',' 'NR==2 {print $2}' $r); echo "${r%%-*} ${f1}" >> all_mmgbsa_scores_pt1.mm; done; sort -k2 -n all_mmgbsa_scores_pt1.mm >> sorted_all_mmgbsa_scores_pt1.mm; head -n 3 sorted_all_mmgbsa_scores_pt1.mm >> selected2_mmgbsa_pt1.mm;

for e3 in e3_library/*pdb; do e=${e3##*/}; while read line; do id=$(echo "$line" | cut -d ' ' -f1 | rev | cut -d '_' -f5-6 | rev); mkdir ${e%.*}_${id}_pair; cp $e3 ${e%.*}_${id}_pair;cp out_nhyd99_mae_amb_${poi%.*}_${id}_linker_docked_short.pdb ${e%.*}_${id}_pair; cp utilities/ProPOSE* ${e%.*}_${id}_pair; cp utilities/*sh ${e%.*}_${id}_pair; cp utilities/*py ${e%.*}_${id}_pair; echo ${poi%.*}_${id}_linker_docked_short.pdb;done <selected1_sorted_filtered_all_rmsd.rmsdum; done;

rm selected1_sorted_filtered_all_rmsd.rmsdum sorted_all_rmsd.rmsdum all_rmsd.rmsdum
for i in *rmsdum; do mv ${i} ${i%.*}.rms;done

for p in *_pair; do (cd $p && sbatch bash_propose.sh *_linker_docked_short.pdb e3*pdb);done

echo "DONE ProPOSE docking; wait for the rest of the pipeline to complete: monitor on the slurm"

mkdir temp_process top_20_pooled process process_2 process_3 process_4 stg2_process ml_scores stg4

echo "temp_process top_20_pooled process process_2 process_3 process_4 stg2_process ml_scores" | xargs -n 1 cp -v ./utilities/{*py,*sh,*csv}

#rm -r selected1_sorted_filtered_all_rmsd.rms selected1_sorted_all_rmsd.rms sorted_all_rmsd.rms all_mmgbsa_scores_pt1.mm sorted_all_mmgbsa_scores_pt1.mm selected2_mmgbsa_pt1.mm *_mmgbsa_rescored.pdb *-out-out.maegz *.mae *out.csv


while [[ $(squeue --name PP_dock | wc -l) -gt 1 ]]; do echo "waiting for protein-protein dock"; sleep 10m; done;
while [[ $(squeue --name LEAP | wc -l) -gt 1 ]]; do echo "proliferating docked poses..."; sleep 10m; done;
while [[ $(squeue --name stage1 | wc -l) -gt 1 ]]; do echo "2nd stage minimization ..."; sleep 10m; done;
while [[ $(squeue --name stage2 | wc -l) -gt 1 ]]; do echo "3rd stage minimization ... "; sleep 10m; done;
while [[ $(squeue --name stage3 | wc -l) -gt 1 ]]; do echo "4th stage minimization ..."; sleep 10m; done;
while [[ $(squeue --name stage4 | wc -l) -gt 1 ]]; do echo "5th stage minimization ... "; sleep 10m; done;
while [[ $(squeue --name PP_dock | wc -l) -gt 1 ]] || [[ $(squeue --name LEAP | wc -l) -gt 1 ]] || [[ $(squeue --name stage1 | wc -l) -gt 1 ]] || [[ $(squeue --name stage2 | wc -l) -gt 1 ]] || [[ $(squeue --name stage3 | wc -l) -gt 1 ]] || [[ $(squeue --name stage4 | wc -l) -gt 1 ]] ; do echo "5th stage minimization ... "; sleep 10m; done;


cd ${genesis}
cp utilities/fetch_results_dummy.sh ./

sbatch fetch_results_dummy.sh ${first_iteration}

while [[ $(squeue --name results | wc -l) -gt 1 ]]; do echo "fetching scores for ML prediction ..."; sleep 10m; done;

cd ${genesis}

#out_nhyd99_mae_complex_3mxf_42_6P3D_1_pp_model_xx04.pdb

mkdir chosen_e3
#e3_ciap1_3cbs_2_docked_rec.pdb
for result in ${first_iteration}/ml_scores/results/*pdb; do id=$(echo "$result" | rev | cut -d '/' -f1 | rev | cut -d '_' -f5-6); echo $result; cp utilities/e3_library/*${id}_*rec.pdb chosen_e3;done

cp -r chosen_e3 ${first_iteration}

cd ${first_iteration}
mkdir dummy_results
mv process_4/*out-out*csv dummy_results
rm -rf *pair temp_process top_20_pooled process process_2 process_3 process_4 stg2_process ml_scores stg4

for e3 in chosen_e3/*pdb; do e=${e3##*/}; while read line; do id=$(echo "$line" | cut -d ' ' -f1 | rev | cut -d '_' -f5-6 | rev); mkdir ${e%.*}_${id}_pair; cp $e3 ${e%.*}_${id}_pair;cp out_nhyd99_mae_amb_${poi%.*}_${id}_linker_docked_short.pdb ${e%.*}_${id}_pair; cp utilities/ProPOSE* ${e%.*}_${id}_pair; cp utilities/*sh ${e%.*}_${id}_pair; cp utilities/*py ${e%.*}_${id}_pair; echo ${poi%.*}_${id}_linker_docked_short.pdb;done <selected1_sorted_filtered_all_rmsd.rms; done;

for p in *_pair; do (cd $p && sbatch bash_propose.sh *_linker_docked_short.pdb e3*pdb);done

echo "DONE ProPOSE docking; wait for the rest of the pipeline to complete: monitor on the slurm"

mkdir temp_process top_20_pooled process process_2 process_3 process_4 stg2_process ml_scores stg4

echo "temp_process top_20_pooled process process_2 process_3 process_4 stg2_process ml_scores" | xargs -n 1 cp -v ./utilities/{*py,*sh,*csv}

#rm -r selected1_sorted_filtered_all_rmsd.rms selected1_sorted_all_rmsd.rms sorted_all_rmsd.rms all_mmgbsa_scores_pt1.mm sorted_all_mmgbsa_scores_pt1.mm selected2_mmgbsa_pt1.mm *_mmgbsa_rescored.pdb *-out-out.maegz *.mae *out.csv


while [[ $(squeue --name PP_dock | wc -l) -gt 1 ]]; do echo "waiting for protein-protein dock"; sleep 10m; done;
while [[ $(squeue --name LEAP | wc -l) -gt 1 ]]; do echo "proliferating docked poses..."; sleep 10m; done;
while [[ $(squeue --name stage1 | wc -l) -gt 1 ]]; do echo "2nd stage minimization ..."; sleep 10m; done;
while [[ $(squeue --name stage2 | wc -l) -gt 1 ]]; do echo "3rd stage minimization ... "; sleep 10m; done;
while [[ $(squeue --name stage3 | wc -l) -gt 1 ]]; do echo "4th stage minimization ..."; sleep 10m; done;
while [[ $(squeue --name stage4 | wc -l) -gt 1 ]]; do echo "5th stage minimization ... "; sleep 10m; done;
while [[ $(squeue --name PP_dock | wc -l) -gt 1 ]] || [[ $(squeue --name LEAP | wc -l) -gt 1 ]] || [[ $(squeue --name stage1 | wc -l) -gt 1 ]] || [[ $(squeue --name stage2 | wc -l) -gt 1 ]] || [[ $(squeue --name stage3 | wc -l) -gt 1 ]] || [[ $(squeue --name stage4 | wc -l) -gt 1 ]] ; do echo "5th stage minimization ... "; sleep 10m; done;


cd ${genesis}

cp utilities/fetch_results.sh ./

sbatch fetch_results.sh ${first_iteration}
while [[ $(squeue --name results | wc -l) -gt 1 ]]; do echo "fetching scores for ML prediction ..."; sleep 10m; done;

cd ${genesis}

#ligandPoi=$1 #poi=$2 #unique_idx=$3 #fragt=$4 #fragt_idx=$5 #fragth=$6 #fragth_idx=$7

### SECOND ITERATION ###

mkdir ${4%.*}_${2%.*}_${5%.*}

cd ${4%.*}_${2%.*}_${5%.*}

second_iteration="${genesis}/${4%.*}_${2%.*}_${5%.*}"

echo "genesis directory is ${genesis}"
echo "2nd_iteration directory is ${second_iteration}"

poi=$2
rm -r -f *rms *mm *log *mol2 *sdf *maps* sed* *cp* linker_library/temp/
cp -r ../utilities/* ./
cp -r ../utilities/ ./

cp -r ${first_iteration}/linker_library/ ./
cp -r ${first_iteration}/conect_cp/ ./
cp -r ../chosen_e3 ./
cp ../{$4,$2} ./

cp $poi ${poi%.*}_cp.pdb
pdb="${poi%.*}_cp.pdb"

## 1) dock ligand to POI

pdb_base=${pdb%.*};
grep "^HETATM" $pdb > ${pdb_base}_ligand.pdb;
ligand=$(awk 'NR==1 {print $4}' ${pdb_base}_ligand.pdb | awk '{print $0}');
sed -i "/HETATM/s/${ligand}/UNL/" $pdb
rm ${pdb_base}_ligand.pdb;
lett="A B C D E F G H"; for L in $lett; do echo $L; sed -i "/ UNL /{s/ UNL ${L} / UNL   /}" $pdb;done

$ICMHOME/icm -vlscluster $ICMHOME/_dockBatch $pdb ligmol=unl
mkdir ${pdb_base}_maps/;
mv -t ${pdb_base}_maps/ *map *htm *dtb *ob;

echo "initiating docking for .." ${fragt};
up_pdb=${pdb_base^^};

sbatch --wait icm_poi.sh ${pdb_base}_maps/D_${up_pdb} ${fragt} ${fragt%.*}_docked;

for i in ${fragt%.*}_docked*sdf; do mv $i ${fragt%.*}_docked.sdf;done

## 2) add ligand to receptor

# remove existing ligand
grep "^HETATM" $pdb > ${pdb_base}_ligand.pdb;
#ligand=$(awk 'NR==1 {print $4}' ${pdb_base}_ligand.pdb | awk '{print tolower($0)}');
rm ${pdb_base}_ligand.pdb;
ligand="UNL"
up_ligand=${ligand^^};
sed -i "/HETATM/{/${up_ligand}/d}" $pdb;
sed -i "/HETATM/{/${ligand}/d}" $pdb;


# add new ligand
cp $pdb ${pdb_base}_rec.pdb;
sed -i '/CONECT/d' ${pdb_base}_rec.pdb
sed -i '/END/d' ${pdb_base}_rec.pdb
obabel -isd ${fragt%.*}_docked.sdf -opdb -O ${fragt%.*}_docked.pdb
cat ${fragt%.*}_docked.pdb >> ${pdb_base}_rec.pdb;
sed -i '/AUTHOR/d' ${pdb_base}_rec.pdb
sed -i '/COMPND/d' ${pdb_base}_rec.pdb

## 3) dock linker on top of the ligand

obabel -ipdb ${fragt%.*}_docked.pdb -omol2 -O ${fragt%.*}_docked.mol2

mol2_lig=${fragt}_docked.mol2
rm -r *_maps
$ICMHOME/icm -vlscluster $ICMHOME/_dockBatch ${pdb_base}_rec.pdb ligmol=${mol2_lig};
mkdir ${pdb_base}_maps/;
mv -t ${pdb_base}_maps/ *map *htm *dtb *ob;

## 4) modify conect records for new POI

cp -r conect/ conect_cp/

for file in conect_cp/*idx; do
        line=$(head -n 1 $file)
        C1=$(cut -d',' -f1 <<<$line)
        C2=$(cut -d',' -f2 <<<$line)
        C3=$(cut -d',' -f3 <<<$line)
        C4=$(cut -d',' -f4 <<<$line)
        C5=$(cut -d',' -f5 <<<$line)
        C6=$(cut -d',' -f6 <<<$line)
        original="$C1,$C2,$C3,$C4,$C5,$C6"
        C1=${fragt_idx}
        if [ -z "$C2" ]; then
                C5=${fragt_idx}
        else
                :
        fi
        new="$C1,$C2,$C3,$C4,$C5,$C6"
        sed -i "s/${original}/${new}/" ${file}
done;

poi=$(ls *_cp_rec.pdb)

for linker in linker_library/temp/*short.pdb; do
link=${linker##*/}

v1=$(echo "$link" | cut -d '_' -f 1)
v2=$(echo "$link" | cut -d '_' -f 2)
v3=${v1}_${v2}

con=$(ls conect_cp/${v3}_*idx)

## amending pdb for linker ligation ## LINKER FIRST ##
line=$(head -n 1 $con)
C1=$(cut -d',' -f1 <<<$line)
C2=$(cut -d',' -f2 <<<$line)
C3=$(cut -d',' -f3 <<<$line)
C4=$(cut -d',' -f4 <<<$line)
C5=$(cut -d',' -f5 <<<$line)
C6=$(cut -d',' -f6 <<<$line)


echo "${poi},${linker},${con},${C1},${C2},${C3},${C4},${C5},${C6}"
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
fi


sbatch stage1.sh $poi ${PWD}/${linker}
done;

while [[ $(squeue --name minimize1 | wc -l) -gt 1 ]]; do echo "waiting for 2nd iteration 1st stage minimization for linkers to finish"; sleep 10s; done;


rm -r *mol2 *maps*
echo "ca y est"

poi=$(ls *_cp_rec.pdb)



for r in *rms; do f1=$(head -n 1 $r | cut -d ' ' -f3); if [ ! -z "$f1" ]; then echo "${r} ${f1}" >>all_rmsd.rms; else : ; fi; done;
sed -i '/inf/d' all_rmsd.rms;

sort -k3 -n all_rmsd.rms >> sorted_all_rmsd.rms;
head -n 250 sorted_all_rmsd.rms| shuf | head -n 100 >> selected1_sorted_all_rmsd.rms;
awk '{print $1,$2}' selected1_sorted_all_rmsd.rms >> selected1_sorted_filtered_all_rmsd.rms;

for e3 in chosen_e3/*pdb; do e=${e3##*/}; while read line; do id=$(echo "$line" | cut -d ' ' -f1 | rev | cut -d '_' -f5-6 | rev); mkdir ${e%.*}_${id}_pair; cp $e3 ${e%.*}_${id}_pair;cp out_nhyd99_mae_amb_${poi%.*}_${id}_linker_docked_short.pdb ${e%.*}_${id}_pair; cp utilities/ProPOSE* ${e%.*}_${id}_pair; cp utilities/*sh ${e%.*}_${id}_pair; cp utilities/*py ${e%.*}_${id}_pair; echo ${poi%.*}_${id}_linker_docked_short.pdb;done <selected1_sorted_filtered_all_rmsd.rms; done;

for p in *_pair; do (cd $p && sbatch bash_propose.sh *_linker_docked_short.pdb e3*pdb);done

echo "DONE ProPOSE docking; wait for the rest of the pipeline to complete: monitor on the slurm"

mkdir temp_process top_20_pooled process process_2 process_3 process_4 stg2_process ml_scores stg4

echo "temp_process top_20_pooled process process_2 process_3 process_4 stg2_process ml_scores" | xargs -n 1 cp -v ./utilities/{*py,*sh,*csv}

#rm -r selected1_sorted_filtered_all_rmsd.rms selected1_sorted_all_rmsd.rms sorted_all_rmsd.rms all_mmgbsa_scores_pt1.mm sorted_all_mmgbsa_scores_pt1.mm selected2_mmgbsa_pt1.mm *_mmgbsa_rescored.pdb *-out-out.maegz *.mae *out.csv


while [[ $(squeue --name PP_dock | wc -l) -gt 1 ]]; do echo "waiting for protein-protein dock"; sleep 10m; done;
while [[ $(squeue --name LEAP | wc -l) -gt 1 ]]; do echo "proliferating docked poses..."; sleep 10m; done;
while [[ $(squeue --name stage1 | wc -l) -gt 1 ]]; do echo "2nd stage minimization ..."; sleep 10m; done;
while [[ $(squeue --name stage2 | wc -l) -gt 1 ]]; do echo "3rd stage minimization ... "; sleep 10m; done;
while [[ $(squeue --name stage3 | wc -l) -gt 1 ]]; do echo "4th stage minimization ..."; sleep 10m; done;
while [[ $(squeue --name stage4 | wc -l) -gt 1 ]]; do echo "5th stage minimization ... "; sleep 10m; done;
while [[ $(squeue --name PP_dock | wc -l) -gt 1 ]] || [[ $(squeue --name LEAP | wc -l) -gt 1 ]] || [[ $(squeue --name stage1 | wc -l) -gt 1 ]] || [[ $(squeue --name stage2 | wc -l) -gt 1 ]] || [[ $(squeue --name stage3 | wc -l) -gt 1 ]] || [[ $(squeue --name stage4 | wc -l) -gt 1 ]] ; do echo "5th stage minimization ... "; sleep 10m; done;


cd ${genesis}
cp utilities/fetch_results.sh ./
sbatch fetch_results.sh ${second_iteration}

while [[ $(squeue --name results | wc -l) -gt 1 ]]; do echo "fetching scores for ML prediction ..."; sleep 10m; done;

cd ${genesis}

### THIRD ITERATION ###

mkdir ${6%.*}_${2%.*}_${7%.*}

cd ${6%.*}_${2%.*}_${7%.*}

third_iteration="${genesis}/${6%.*}_${2%.*}_${7%.*}"

echo "genesis directory is ${genesis}"
echo "3rd_iteration directory is ${third_iteration}"

poi=$2
rm -r -f *rms *mm *log *mol2 *sdf *maps* sed* *cp* linker_library/temp/
cp -r ../utilities/* ./
cp -r ../utilities/ ./

cp -r ${first_iteration}/linker_library/ ./
cp -r ${first_iteration}/conect_cp/ ./
cp -r ../chosen_e3 ./
cp ../{$6,$2} ./

cp $poi ${poi%.*}_cp.pdb
pdb="${poi%.*}_cp.pdb"

## 1) dock ligand to POI

pdb_base=${pdb%.*};
grep "^HETATM" $pdb > ${pdb_base}_ligand.pdb;
ligand=$(awk 'NR==1 {print $4}' ${pdb_base}_ligand.pdb | awk '{print $0}');
sed -i "/HETATM/s/${ligand}/UNL/" $pdb
rm ${pdb_base}_ligand.pdb;
lett="A B C D E F G H"; for L in $lett; do echo $L; sed -i "/ UNL /{s/ UNL ${L} / UNL   /}" $pdb;done

$ICMHOME/icm -vlscluster $ICMHOME/_dockBatch $pdb ligmol=unl
mkdir ${pdb_base}_maps/;
mv -t ${pdb_base}_maps/ *map *htm *dtb *ob;

echo "initiating docking for .." $fragth;
up_pdb=${pdb_base^^};

sbatch --wait icm_poi.sh ${pdb_base}_maps/D_${up_pdb} $fragth ${fragth%.*}_docked;

for i in ${fragth%.*}_docked*sdf; do mv $i ${fragth%.*}_docked.sdf;done

## 2) add ligand to receptor

# remove existing ligand
grep "^HETATM" $pdb > ${pdb_base}_ligand.pdb;
#ligand=$(awk 'NR==1 {print $4}' ${pdb_base}_ligand.pdb | awk '{print tolower($0)}');
rm ${pdb_base}_ligand.pdb;
ligand="UNL"

up_ligand=${ligand^^};
sed -i "/HETATM/{/${up_ligand}/d}" $pdb;
sed -i "/HETATM/{/${ligand}/d}" $pdb;


# add new ligand
cp $pdb ${pdb_base}_rec.pdb;
sed -i '/CONECT/d' ${pdb_base}_rec.pdb
sed -i '/END/d' ${pdb_base}_rec.pdb
obabel -isd ${fragth%.*}_docked.sdf -opdb -O ${fragth%.*}_docked.pdb
cat ${fragth%.*}_docked.pdb >> ${pdb_base}_rec.pdb;
sed -i '/AUTHOR/d' ${pdb_base}_rec.pdb
sed -i '/COMPND/d' ${pdb_base}_rec.pdb

## 3) dock linker on top of the ligand

obabel -ipdb ${fragth%.*}_docked.pdb -omol2 -O ${fragth%.*}_docked.mol2

mol2_lig=${fragth}_docked.mol2
rm -r *_maps
$ICMHOME/icm -vlscluster $ICMHOME/_dockBatch ${pdb_base}_rec.pdb ligmol=${mol2_lig};
mkdir ${pdb_base}_maps/;
mv -t ${pdb_base}_maps/ *map *htm *dtb *ob;

## 4) modify conect records for new POI

cp -r conect/ conect_cp/

for file in conect_cp/*idx; do
        line=$(head -n 1 $file)
        C1=$(cut -d',' -f1 <<<$line)
        C2=$(cut -d',' -f2 <<<$line)
        C3=$(cut -d',' -f3 <<<$line)
        C4=$(cut -d',' -f4 <<<$line)
        C5=$(cut -d',' -f5 <<<$line)
        C6=$(cut -d',' -f6 <<<$line)
        original="$C1,$C2,$C3,$C4,$C5,$C6"
        C1=${fragth_idx}
        if [ -z "$C2" ]; then
                C5=${fragth_idx}
        else
                :
        fi
        new="$C1,$C2,$C3,$C4,$C5,$C6"
        sed -i "s/${original}/${new}/" ${file}
done;

poi=$(ls *_cp_rec.pdb)

for linker in linker_library/temp/*short.pdb; do
link=${linker##*/}

v1=$(echo "$link" | cut -d '_' -f 1)
v2=$(echo "$link" | cut -d '_' -f 2)
v3=${v1}_${v2}

con=$(ls conect_cp/${v3}_*idx)

## amending pdb for linker ligation ## LINKER FIRST ##
line=$(head -n 1 $con)
C1=$(cut -d',' -f1 <<<$line)
C2=$(cut -d',' -f2 <<<$line)
C3=$(cut -d',' -f3 <<<$line)
C4=$(cut -d',' -f4 <<<$line)
C5=$(cut -d',' -f5 <<<$line)
C6=$(cut -d',' -f6 <<<$line)


echo "${poi},${linker},${con},${C1},${C2},${C3},${C4},${C5},${C6}"
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
fi

sbatch stage1.sh $poi ${PWD}/${linker}
done;

while [[ $(squeue --name minimize1 | wc -l) -gt 1 ]]; do echo "waiting for 2nd iteration 1st stage minimization for linkers to finish"; sleep 10s; done;


rm -r *mol2 *maps*
echo "ca y est"

poi=$(ls *_cp_rec.pdb)



for r in *rms; do f1=$(head -n 1 $r | cut -d ' ' -f3); if [ ! -z "$f1" ]; then echo "${r} ${f1}" >>all_rmsd.rms; else : ; fi; done;
sed -i '/inf/d' all_rmsd.rms;

sort -k3 -n all_rmsd.rms >> sorted_all_rmsd.rms;
head -n 250 sorted_all_rmsd.rms| shuf | head -n 100 >> selected1_sorted_all_rmsd.rms;
awk '{print $1,$2}' selected1_sorted_all_rmsd.rms >> selected1_sorted_filtered_all_rmsd.rms;

for e3 in chosen_e3/*pdb; do e=${e3##*/}; while read line; do id=$(echo "$line" | cut -d ' ' -f1 | rev | cut -d '_' -f5-6 | rev); mkdir ${e%.*}_${id}_pair; cp $e3 ${e%.*}_${id}_pair;cp out_nhyd99_mae_amb_${poi%.*}_${id}_linker_docked_short.pdb ${e%.*}_${id}_pair; cp utilities/ProPOSE* ${e%.*}_${id}_pair; cp utilities/*sh ${e%.*}_${id}_pair; cp utilities/*py ${e%.*}_${id}_pair; echo ${poi%.*}_${id}_linker_docked_short.pdb;done <selected1_sorted_filtered_all_rmsd.rms; done;

for p in *_pair; do (cd $p && sbatch bash_propose.sh *_linker_docked_short.pdb e3*pdb);done

echo "DONE ProPOSE docking; wait for the rest of the pipeline to complete: monitor on the slurm"

mkdir temp_process top_20_pooled process process_2 process_3 process_4 stg2_process ml_scores stg4

echo "temp_process top_20_pooled process process_2 process_3 process_4 stg2_process ml_scores" | xargs -n 1 cp -v ./utilities/{*py,*sh,*csv}

#rm -r selected1_sorted_filtered_all_rmsd.rms selected1_sorted_all_rmsd.rms sorted_all_rmsd.rms all_mmgbsa_scores_pt1.mm sorted_all_mmgbsa_scores_pt1.mm selected2_mmgbsa_pt1.mm *_mmgbsa_rescored.pdb *-out-out.maegz *.mae *out.csv


while [[ $(squeue --name PP_dock | wc -l) -gt 1 ]]; do echo "waiting for protein-protein dock"; sleep 10m; done;
while [[ $(squeue --name LEAP | wc -l) -gt 1 ]]; do echo "proliferating docked poses..."; sleep 10m; done;
while [[ $(squeue --name stage1 | wc -l) -gt 1 ]]; do echo "2nd stage minimization ..."; sleep 10m; done;
while [[ $(squeue --name stage2 | wc -l) -gt 1 ]]; do echo "3rd stage minimization ... "; sleep 10m; done;
while [[ $(squeue --name stage3 | wc -l) -gt 1 ]]; do echo "4th stage minimization ..."; sleep 10m; done;
while [[ $(squeue --name stage4 | wc -l) -gt 1 ]]; do echo "5th stage minimization ... "; sleep 10m; done;
while [[ $(squeue --name PP_dock | wc -l) -gt 1 ]] || [[ $(squeue --name LEAP | wc -l) -gt 1 ]] || [[ $(squeue --name stage1 | wc -l) -gt 1 ]] || [[ $(squeue --name stage2 | wc -l) -gt 1 ]] || [[ $(squeue --name stage3 | wc -l) -gt 1 ]] || [[ $(squeue --name stage4 | wc -l) -gt 1 ]] ; do echo "5th stage minimization ... "; sleep 10m; done;


cd ${genesis}
cp utilities/fetch_results.sh ./
sbatch fetch_results.sh ${third_iteration}
while [[ $(squeue --name results | wc -l) -gt 1 ]]; do echo "fetching scores for ML prediction ..."; sleep 10m; done;

cd ${genesis}

mkdir results

cp */ml_scores/results results


