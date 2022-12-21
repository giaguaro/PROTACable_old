#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=PP_dock

source ~/.bashrc
conda activate amber

echo "%{j}.out"
p=$1
e3=$2

#mkdir ${poi%.*}_${e3%.*}_pp_docking_trial2
#out_nhyd99_mae_amb_poi_mpro_cp_rec_2GQG_35_linker_docked_short.pdb
id=$(echo "$p" | cut -d '.' -f1 | rev | cut -d '_' -f4-5 | rev)

cp $p poi_target_${id}.pdb
poi="poi_target_${id}.pdb"

#mkdir ${poi%.*}_${e3%.*}_pp_docking_trial2

#cp $poi ${poi%.*}_${e3%.*}_pp_docking_trial2
#cp $e3 ${poi%.*}_${e3%.*}_pp_docking_trial2

#cp ProPOSE ${poi%.*}_${e3%.*}_pp_docking_trial2
#cp ProPOSE.lic ${poi%.*}_${e3%.*}_pp_docking_trial2
#cp list_cdr_atoms.sh ${poi%.*}_${e3%.*}_pp_docking_trial2


export PROPOSEHOME=$HOME/propose

#cp $poi ${poi%.*}_${e3%.*}_pp_docking
#cp $e3 ${poi%.*}_${e3%.*}_pp_docking

#cp ProPOSE ${poi%.*}_${e3%.*}_pp_docking
#cp ProPOSE.lic ${poi%.*}_${e3%.*}_pp_docking
#cp list_cdr_atoms.sh ${poi%.*}_${e3%.*}_pp_docking 

#cd ${poi%.*}_${e3%.*}_pp_docking_trial2

#sed -i -e '/HETATM/,/CONECT/!b' -e '/TER/d' $poi
#sed -i -e '/HETATM/,/CONECT/!b' -e '/TER/d' $e3
#sed -i -e '/ATOM/,/CONECT/!b' -e '/TER/d' $poi
#sed -i -e '/ATOM/,/CONECT/!b' -e '/TER/d' $e3


count_ter=$( grep -c TER $poi)
count_ter_n=${count_ter}-1
#if [ "$count_ter" -gt 1 ]; then
#	sed -i '0,/TER/{/TER/d}' $poi
#else
#	:
#fi
#
count_ter2=$( grep -c TER $e3)
count_ter2_n=${count_ter2}-1
#
#if [ "$count_ter2" -gt 1 ]; then
#        sed -i '0,/TER/{/TER/d}' $e3
#else
#        :
#fi

echo "you have ${count_ter_n} break(s) in your chain ${poi%.*}"
echo "you have ${count_ter2_n} break(s) in your chain ${e3%.*}"
#sed -i -e '/TER/,/END/!b' -e 's/ATOM  /HETATM/g' $poi
#sed -i -e '/TER/,/END/!b' -e 's/ATOM  /HETATM/g' $e3

rm -f *scores.csv *special* *hit *trial* *prep* *ligand* sqm.* ANTECHAMBER* ATOMTYPE* *log *noHET* *sed* *mol2 *sslink* *txt

cp $poi ${poi%.*}_prep1.pdb
cp $e3 ${e3%.*}_prep1.pdb


echo "renaming HETATMs"

o_atoms="C O N F S P"
t_atoms="Cl Br"
o_Hatoms="H"
t_Hatoms="HN HO HS HP"
th_Hatoms="HXT HN1 HN2 HN3 HO1 1HO HO2 2HO HO3 3HO HS1 1HS HS2 2SH HS3 3SH HP1 1HP HP2 2HP HP3 3HP"
f_Hatoms="1HXT 2HXT 3HXT"

#POI
####
for i in {1..9}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${poi%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${poi%.*}_prep1.pdb;done;done

for i in {1..9}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${poi%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${poi%.*}_prep1.pdb;done;done

sed -i 's/ H  1 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb; sed -i 's/ H  2 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb; sed -i 's/ H  3 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb;

for hyd in $t_Hatoms; do sed -i "s/ $hyd  UNL/ H   UNL/g" ${poi%.*}_prep1.pdb;done
for hyd in $th_Hatoms; do sed -i "s/ $hyd UNL/ H   UNL/g" ${poi%.*}_prep1.pdb;done
for hyd in $f_Hatoms; do sed -i "s/ $hyd UNL/  H   UNL/g" ${poi%.*}_prep1.pdb;done


for i in {1..9}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i}  UNL/}" ${poi%.*}_prep1.pdb;done
for i in {10..80}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i} UNL/}" ${poi%.*}_prep1.pdb;done

#E3
###
for i in {1..9}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${e3%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${e3%.*}_prep1.pdb;done;done

for i in {1..9}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${e3%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${e3%.*}_prep1.pdb;done;done

sed -i 's/ H  1 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb; sed -i 's/ H  2 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb; sed -i 's/ H  3 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb;

for hyd in $t_Hatoms; do sed -i "s/ $hyd  UNL/ H   UNL/g" ${e3%.*}_prep1.pdb;done
for hyd in $th_Hatoms; do sed -i "s/ $hyd UNL/ H   UNL/g" ${e3%.*}_prep1.pdb;done
for hyd in $f_Hatoms; do sed -i "s/ $hyd UNL/  H   UNL/g" ${e3%.*}_prep1.pdb;done

for i in {1..9}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i}  UNL/}" ${e3%.*}_prep1.pdb;done
for i in {10..80}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i} UNL/}" ${e3%.*}_prep1.pdb;done
####################################################################################################################################################

#### spesh case
cp ${poi%.*}_prep1.pdb ${poi%.*}_prep1_special.pdb
grep "^HETATM" ${poi%.*}_prep1_special.pdb > ${poi%.*}_prep1_special_lig.pdb
grep "^CONECT" ${poi%.*}_prep1_special.pdb >> ${poi%.*}_prep1_special_lig.pdb
grep "^END" ${poi%.*}_prep1_special.pdb >> ${poi%.*}_prep1_special_lig.pdb

sed -i '/HETATM/d' ${poi%.*}_prep1_special.pdb
sed -i '/CONECT/d' ${poi%.*}_prep1_special.pdb
sed -i '/END/d' ${poi%.*}_prep1_special.pdb


pdb4amber -i ${poi%.*}_prep1_special.pdb -o ${poi%.*}_r_special.pdb --nohyd
sed -i '/HETATM/d' ${poi%.*}_r_special.pdb
sed -i '/CONECT/d' ${poi%.*}_r_special.pdb
sed -i '/END/d' ${poi%.*}_r_special.pdb

cat ${poi%.*}_prep1_special_lig.pdb >> ${poi%.*}_r_special.pdb

pdb4amber -i ${poi%.*}_r_special.pdb -o ${poi%.*}_r_special2.pdb
#####

pdb4amber -i ${poi%.*}_prep1.pdb -o ${poi%.*}_r_prep1.pdb

if [ ! -f ${poi%.*}_r_prep1.pdb ] || [ ! -s ${poi%.*}_r_prep1.pdb ]; then
	conda deactivate amber
	conda activate py39
	grep "^HETATM" ${poi%.*}_prep1.pdb > ${poi%.*}_ligand0.pdb
	grep "^CONECT" ${poi%.*}_prep1.pdb >> ${poi%.*}_ligand0.pdb
	grep "^END" ${poi%.*}_prep1.pdb >> ${poi%.*}_ligand0.pdb
	
	sed -i '/HETATM/d' ${poi%.*}_prep1.pdb
	sed -i '/CONECT/d' ${poi%.*}_prep1.pdb
	sed -i '/END/d' ${poi%.*}_prep1.pdb

	pdbfixer ${poi%.*}_prep1.pdb --replace-nonstandard --add-residues --output="fixed_${poi%.*}_prep1.pdb"
	sed -i '/END/d' fixed_${poi%.*}_prep1.pdb
	sed -i '/TER.*UNL\|UNL.*TER/d' fixed_${poi%.*}_prep1.pdb
	cat ${poi%.*}_ligand0.pdb >> fixed_${poi%.*}_prep1.pdb
	
	conda deactivate py39
	conda activate amber

	pdb4amber -i fixed_${poi%.*}_prep1.pdb -o ${poi%.*}_r_prep1.pdb
	cp $poi ${poi%.*}_prep1.pdb
	
	echo "renaming HETATMs"

        o_atoms="C O N F S P"
        t_atoms="Cl Br"
        o_Hatoms="H"
        t_Hatoms="HN HO HS HP"
        th_Hatoms="HXT HN1 HN2 HN3 HO1 1HO HO2 2HO HO3 3HO HS1 1HS HS2 2SH HS3 3SH HP1 1HP HP2 2HP HP3 3HP"
        f_Hatoms="1HXT 2HXT 3HXT"

        #POI
        ####
        for i in {1..9}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${poi%.*}_prep1.pdb;done;done
        for i in {10..80}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${poi%.*}_prep1.pdb;done;done

        for i in {1..9}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${poi%.*}_prep1.pdb;done;done
        for i in {10..80}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${poi%.*}_prep1.pdb;done;done

        sed -i 's/ H  1 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb; sed -i 's/ H  2 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb; sed -i 's/ H  3 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb;

        for hyd in $t_Hatoms; do sed -i "s/ $hyd  UNL/ H   UNL/g" ${poi%.*}_prep1.pdb;done
        for hyd in $th_Hatoms; do sed -i "s/ $hyd UNL/ H   UNL/g" ${poi%.*}_prep1.pdb;done
        for hyd in $f_Hatoms; do sed -i "s/ $hyd UNL/  H   UNL/g" ${poi%.*}_prep1.pdb;done


        for i in {1..9}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i}  UNL/}" ${poi%.*}_prep1.pdb;done
        for i in {10..80}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i} UNL/}" ${poi%.*}_prep1.pdb;done


	cp fixed_${poi%.*}_prep1.pdb ${poi%.*}_prep1_special.pdb
	## spesh case
	rm ${poi%.*}_r_special.pdb ${poi%.*}_r_special2.pdb 

	grep "^HETATM" ${poi%.*}_prep1_special.pdb > ${poi%.*}_prep1_special_lig.pdb
	grep "^CONECT" ${poi%.*}_prep1_special.pdb >> ${poi%.*}_prep1_special_lig.pdb
	grep "^END" ${poi%.*}_prep1_special.pdb >> ${poi%.*}_prep1_special_lig.pdb

	sed -i '/HETATM/d' ${poi%.*}_prep1_special.pdb
	sed -i '/CONECT/d' ${poi%.*}_prep1_special.pdb
	sed -i '/END/d' ${poi%.*}_prep1_special.pdb


	pdb4amber -i ${poi%.*}_prep1_special.pdb -o ${poi%.*}_r_special.pdb --nohyd
	sed -i '/HETATM/d' ${poi%.*}_r_special.pdb
	sed -i '/CONECT/d' ${poi%.*}_r_special.pdb
	sed -i '/END/d' ${poi%.*}_r_special.pdb

	cat ${poi%.*}_prep1_special_lig.pdb >> ${poi%.*}_r_special.pdb

	pdb4amber -i ${poi%.*}_r_special.pdb -o ${poi%.*}_r_special2.pdb
fi

pdb4amber -i ${e3%.*}_prep1.pdb -o ${e3%.*}_r_prep1.pdb



#### spesh case
cp ${e3%.*}_prep1.pdb ${e3%.*}_prep1_special.pdb
grep "^HETATM" ${e3%.*}_prep1_special.pdb > ${e3%.*}_prep1_special_lig.pdb
grep "^CONECT" ${e3%.*}_prep1_special.pdb >> ${e3%.*}_prep1_special_lig.pdb
grep "^END" ${e3%.*}_prep1_special.pdb >> ${e3%.*}_prep1_special_lig.pdb

sed -i '/HETATM/d' ${e3%.*}_prep1_special.pdb
sed -i '/CONECT/d' ${e3%.*}_prep1_special.pdb
sed -i '/END/d' ${e3%.*}_prep1_special.pdb


pdb4amber -i ${e3%.*}_prep1_special.pdb -o ${e3%.*}_r_special.pdb --nohyd
sed -i '/HETATM/d' ${e3%.*}_r_special.pdb
sed -i '/CONECT/d' ${e3%.*}_r_special.pdb
sed -i '/END/d' ${e3%.*}_r_special.pdb

cat ${e3%.*}_prep1_special_lig.pdb >> ${e3%.*}_r_special.pdb

pdb4amber -i ${e3%.*}_r_special.pdb -o ${e3%.*}_r_special2.pdb
#####

if [ ! -f ${e3%.*}_r_prep1.pdb ] || [ ! -s ${e3%.*}_r_prep1.pdb ]; then
        conda deactivate amber
        conda activate py39
        grep "^HETATM" ${e3%.*}_prep1.pdb > ${e3%.*}_ligand0.pdb
        grep "^CONECT" ${e3%.*}_prep1.pdb >> ${e3%.*}_ligand0.pdb
        grep "^END" ${e3%.*}_prep1.pdb >> ${e3%.*}_ligand0.pdb

        sed -i '/HETATM/d' ${e3%.*}_prep1.pdb
        sed -i '/CONECT/d' ${e3%.*}_prep1.pdb
        sed -i '/END/d' ${e3%.*}_prep1.pdb

        pdbfixer ${e3%.*}_prep1.pdb --replace-nonstandard --add-residues --output="fixed_${e3%.*}_prep1.pdb"
        sed -i '/END/d' fixed_${e3%.*}_prep1.pdb
        sed -i '/TER.*UNL\|UNL.*TER/d' fixed_${e3%.*}_prep1.pdb
        cat ${e3%.*}_ligand0.pdb >> fixed_${e3%.*}_prep1.pdb

        conda deactivate py39
        conda activate amber

        pdb4amber -i fixed_${e3%.*}_prep1.pdb -o ${e3%.*}_r_prep1.pdb
        cp $e3 ${e3%.*}_prep1.pdb
	

        #E3
        ###
        for i in {1..9}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${e3%.*}_prep1.pdb;done;done
        for i in {10..80}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${e3%.*}_prep1.pdb;done;done

        for i in {1..9}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${e3%.*}_prep1.pdb;done;done
        for i in {10..80}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${e3%.*}_prep1.pdb;done;done

        sed -i 's/ H  1 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb; sed -i 's/ H  2 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb; sed -i 's/ H  3 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb;

        for hyd in $t_Hatoms; do sed -i "s/ $hyd  UNL/ H   UNL/g" ${e3%.*}_prep1.pdb;done
        for hyd in $th_Hatoms; do sed -i "s/ $hyd UNL/ H   UNL/g" ${e3%.*}_prep1.pdb;done
        for hyd in $f_Hatoms; do sed -i "s/ $hyd UNL/  H   UNL/g" ${e3%.*}_prep1.pdb;done

        for i in {1..9}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i}  UNL/}" ${e3%.*}_prep1.pdb;done
        for i in {10..80}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i} UNL/}" ${e3%.*}_prep1.pdb;done
	
        cp fixed_${e3%.*}_prep1.pdb ${e3%.*}_prep1_special.pdb
        ## spesh case
        rm ${e3%.*}_r_special.pdb ${e3%.*}_r_special2.pdb

        grep "^HETATM" ${e3%.*}_prep1_special.pdb > ${e3%.*}_prep1_special_lig.pdb
        grep "^CONECT" ${e3%.*}_prep1_special.pdb >> ${e3%.*}_prep1_special_lig.pdb
        grep "^END" ${e3%.*}_prep1_special.pdb >> ${e3%.*}_prep1_special_lig.pdb

        sed -i '/HETATM/d' ${e3%.*}_prep1_special.pdb
        sed -i '/CONECT/d' ${e3%.*}_prep1_special.pdb
        sed -i '/END/d' ${e3%.*}_prep1_special.pdb


        pdb4amber -i ${e3%.*}_prep1_special.pdb -o ${e3%.*}_r_special.pdb --nohyd
        sed -i '/HETATM/d' ${e3%.*}_r_special.pdb
        sed -i '/CONECT/d' ${e3%.*}_r_special.pdb
        sed -i '/END/d' ${e3%.*}_r_special.pdb

        cat ${e3%.*}_prep1_special_lig.pdb >> ${e3%.*}_r_special.pdb

        pdb4amber -i ${e3%.*}_r_special.pdb -o ${e3%.*}_r_special2.pdb
fi


echo "renaming HETATMs"

o_atoms="C O N F S P"
t_atoms="Cl Br"
o_Hatoms="H"
t_Hatoms="HN HO HS HP"
th_Hatoms="HXT HN1 HN2 HN3 HO1 1HO HO2 2HO HO3 3HO HS1 1HS HS2 2SH HS3 3SH HP1 1HP HP2 2HP HP3 3HP"
f_Hatoms="1HXT 2HXT 3HXT"

#POI
####
for i in {1..9}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${poi%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${poi%.*}_prep1.pdb;done;done

for i in {1..9}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${poi%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${poi%.*}_prep1.pdb;done;done

sed -i 's/ H  1 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb; sed -i 's/ H  2 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb; sed -i 's/ H  3 UNL/  H   UNL/g' ${poi%.*}_prep1.pdb;

for hyd in $t_Hatoms; do sed -i "s/ $hyd  UNL/ H   UNL/g" ${poi%.*}_prep1.pdb;done
for hyd in $th_Hatoms; do sed -i "s/ $hyd UNL/ H   UNL/g" ${poi%.*}_prep1.pdb;done
for hyd in $f_Hatoms; do sed -i "s/ $hyd UNL/  H   UNL/g" ${poi%.*}_prep1.pdb;done


for i in {1..9}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i}  UNL/}" ${poi%.*}_prep1.pdb;done
for i in {10..80}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i} UNL/}" ${poi%.*}_prep1.pdb;done

#E3
###
for i in {1..9}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${e3%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $o_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${e3%.*}_prep1.pdb;done;done

for i in {1..9}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i}  UNL/}" ${e3%.*}_prep1.pdb;done;done
for i in {10..80}; do for atm in $t_atoms; do sed -i "$((i-1)),/ $atm   UNL/{s/ $atm   UNL/ ${atm}${i} UNL/}" ${e3%.*}_prep1.pdb;done;done

sed -i 's/ H  1 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb; sed -i 's/ H  2 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb; sed -i 's/ H  3 UNL/  H   UNL/g' ${e3%.*}_prep1.pdb;

for hyd in $t_Hatoms; do sed -i "s/ $hyd  UNL/ H   UNL/g" ${e3%.*}_prep1.pdb;done
for hyd in $th_Hatoms; do sed -i "s/ $hyd UNL/ H   UNL/g" ${e3%.*}_prep1.pdb;done
for hyd in $f_Hatoms; do sed -i "s/ $hyd UNL/  H   UNL/g" ${e3%.*}_prep1.pdb;done

for i in {1..9}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i}  UNL/}" ${e3%.*}_prep1.pdb;done
for i in {10..80}; do sed -i "$((i-1)),/ H   UNL/{s/ H   UNL/ H${i} UNL/}" ${e3%.*}_prep1.pdb;done


#if [ ! -f ${poi%.*}_r_prep1.pdb ] || [ ! -s ${poi%.*}_r_prep1.pdb ]; then
#        obabel -ipdb ${poi%.*}_prep1.pdb -opdb -O ${poi%.*}_r_prep1_ob.pdb
#	sed -i '/CONECT/d' ${poi%.*}_r_prep1_ob.pdb
#	sed -i '/MASTER/d' ${poi%.*}_r_prep1_ob.pdb
#	pdb4amber -i i ${poi%.*}_r_prep1_ob.pdb -o ${poi%.*}_r_prep1.pdb
#else
#        :
#fi

#if [ ! -f ${e3%.*}_r_prep1.pdb ] || [ ! -s ${e3%.*}_r_prep1.pdb ]; then
#        obabel -ipdb ${e3%.*}_prep1.pdb -opdb -O ${e3%.*}_r_prep1_ob.pdb
#        sed -i '/CONECT/d' ${e3%.*}_r_prep1_ob.pdb
#        sed -i '/MASTER/d' ${e3%.*}_r_prep1_ob.pdb
#        pdb4amber -i i ${e3%.*}_r_prep1_ob.pdb -o ${e3%.*}_r_prep1.pdb
#else
#        :
#fi

echo "initiating docking for complexes: $poi and $e3"

grep "^HETATM" ${poi%.*}_prep1.pdb > ${poi%.*}_ligand.pdb
#grep "^CONECT" $poi >> ${poi%.*}_ligand.pdb
grep "^HETATM" ${e3%.*}_prep1.pdb > ${e3%.*}_ligand.pdb
#grep "^CONECT" $e3 >> ${e3%.*}_ligand.pdb

sed -i '/UNL/!d' ${poi%.*}_ligand.pdb

sed -i '/UNL/!d' ${e3%.*}_ligand.pdb


sed -i '/ACE/d' ${poi%.*}_r_prep1.pdb
sed -i '/ACE/d' ${e3%.*}_r_prep1.pdb
sed -i '/NME/d' ${poi%.*}_r_prep1.pdb
sed -i '/NME/d' ${e3%.*}_r_prep1.pdb

sed -i '/HETATM/s/CL/Cl/g' ${poi%.*}_r_prep1.pdb
sed -i '/HETATM/s/CL/Cl/g' ${e3%.*}_r_prep1.pdb
sed -i '/HETATM/s/BR/Br/g' ${poi%.*}_r_prep1.pdb
sed -i '/HETATM/s/BR/Br/g' ${e3%.*}_r_prep1.pdb

############################################################################################################################################################
sed -i 's/CSD/CYS/g' ${poi%.*}_r_prep1.pdb                 #3-SULFINOALANINE
sed -i 's/HYP/PRO/g' ${poi%.*}_r_prep1.pdb                 #4-HYDROXYPROLINE
sed -i 's/BMT/THR/g' ${poi%.*}_r_prep1.pdb                 #4-METHYL-4-[(E)-2-BUTENYL]-4,N-METHYL-THREONINE
sed -i 's/5HP/GLU/g' ${poi%.*}_r_prep1.pdb                 #5-HYDROXYPROLINE
sed -i 's/ABA/ALA/g' ${poi%.*}_r_prep1.pdb                 #ALPHA-AMINOBUTYRIC_ACID
sed -i 's/AIB/ALA/g' ${poi%.*}_r_prep1.pdb                 #ALPHA-AMINOISOBUTYRIC_ACID
sed -i 's/CSW/CYS/g' ${poi%.*}_r_prep1.pdb                 #CYSTEINE-S-DIOXIDE
sed -i 's/OCS/CYS/g' ${poi%.*}_r_prep1.pdb                 #CYSTEINESULFONIC_ACID
sed -i 's/DAL/ALA/g' ${poi%.*}_r_prep1.pdb                 #D-ALANINE
sed -i 's/DAR/ARG/g' ${poi%.*}_r_prep1.pdb                 #D-ARGININE
sed -i 's/DSG/ASN/g' ${poi%.*}_r_prep1.pdb                 #D-ASPARAGINE
sed -i 's/DSP/ASP/g' ${poi%.*}_r_prep1.pdb                 #D-ASPARTATE
sed -i 's/DCY/CYS/g' ${poi%.*}_r_prep1.pdb                 #D-CYSTEINE
sed -i 's/CRO/CRO/g' ${poi%.*}_r_prep1.pdb                 #DECARBOXY(PARAHYDROXYBENZYLIDENE-IMIDAZOLIDINONE)THREONINE
sed -i 's/DGL/GLU/g' ${poi%.*}_r_prep1.pdb                 #D-GLUTAMATE
sed -i 's/DGN/GLN/g' ${poi%.*}_r_prep1.pdb                 #D-GLUTAMINE
sed -i 's/DHI/HIS/g' ${poi%.*}_r_prep1.pdb                 #D-HISTIDINE
sed -i 's/DIL/ILE/g' ${poi%.*}_r_prep1.pdb                 #D-ISOLEUCINE
sed -i 's/DIV/VAL/g' ${poi%.*}_r_prep1.pdb                 #D-ISOVALINE
sed -i 's/DLE/LEU/g' ${poi%.*}_r_prep1.pdb                 #D-LEUCINE
sed -i 's/DLY/LYS/g' ${poi%.*}_r_prep1.pdb                 #D-LYSINE
sed -i 's/DPN/PHE/g' ${poi%.*}_r_prep1.pdb                 #D-PHENYLALANINE
sed -i 's/DPR/PRO/g' ${poi%.*}_r_prep1.pdb                 #D-PROLINE
sed -i 's/DSN/SER/g' ${poi%.*}_r_prep1.pdb                 #D-SERINE
sed -i 's/DTH/THR/g' ${poi%.*}_r_prep1.pdb                 #D-THREONINE
sed -i 's/DTR/DTR/g' ${poi%.*}_r_prep1.pdb                 #D-TRYPTOPHANE
sed -i 's/DTY/TYR/g' ${poi%.*}_r_prep1.pdb                 #D-TYROSINE
sed -i 's/DVA/VAL/g' ${poi%.*}_r_prep1.pdb                 #D-VALINE
sed -i 's/CGU/GLU/g' ${poi%.*}_r_prep1.pdb                 #GAMMA-CARBOXY-GLUTAMIC_ACID
sed -i 's/KCX/LYS/g' ${poi%.*}_r_prep1.pdb                 #LYSINE_NZ-CARBOXYLIC_ACID
sed -i 's/LLP/LYS/g' ${poi%.*}_r_prep1.pdb                 #LYSINE-PYRIDOXAL-5'-PHOSPHATE
sed -i 's/CXM/MET/g' ${poi%.*}_r_prep1.pdb                 #N-CARBOXYMETHIONINE
sed -i 's/FME/MET/g' ${poi%.*}_r_prep1.pdb                 #N-FORMYLMETHIONINE
sed -i 's/MLE/LEU/g' ${poi%.*}_r_prep1.pdb                 #N-METHYLLEUCINE
sed -i 's/MVA/VAL/g' ${poi%.*}_r_prep1.pdb                 #N-METHYLVALINE
sed -i 's/NLE/LEU/g' ${poi%.*}_r_prep1.pdb                 #NORLEUCINE
sed -i 's/PTR/TYR/g' ${poi%.*}_r_prep1.pdb                 #O-PHOSPHOTYROSINE
sed -i 's/ORN/ALA/g' ${poi%.*}_r_prep1.pdb                 #ORNITHINE
sed -i 's/SEP/SER/g' ${poi%.*}_r_prep1.pdb                 #PHOSPHOSERINE
sed -i 's/TPO/THR/g' ${poi%.*}_r_prep1.pdb                 #PHOSPHOTHREONINE
sed -i 's/PCA/GLU/g' ${poi%.*}_r_prep1.pdb                 #PYROGLUTAMIC_ACID
sed -i 's/SAR/GLY/g' ${poi%.*}_r_prep1.pdb                 #SARCOSINE
sed -i 's/CEA/CYS/g' ${poi%.*}_r_prep1.pdb                 #S-HYDROXY-CYSTEINE
sed -i 's/CSO/CYS/g' ${poi%.*}_r_prep1.pdb                 #S-HYDROXYCYSTEINE
sed -i 's/CSS/CYS/g' ${poi%.*}_r_prep1.pdb                 #S-MERCAPTOCYSTEINE
sed -i 's/CSX/CYS/g' ${poi%.*}_r_prep1.pdb                 #S-OXY_CYSTEINE
sed -i 's/CME/CYS/g' ${poi%.*}_r_prep1.pdb                 #S,S-(2-HYDROXYETHYL)THIOCYSTEINE
sed -i 's/TYS/TYR/g' ${poi%.*}_r_prep1.pdb                 #SULFONATED_TYROSINE
sed -i 's/TPQ/PHE/g' ${poi%.*}_r_prep1.pdb                 #TOPO-QUINONE
sed -i 's/STY/TYR/g' ${poi%.*}_r_prep1.pdb                 #TYROSINE-O-SULPHONIC_ACID
sed -i 's/CCS/CYS/g' ${poi%.*}_r_prep1.pdb                 #https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i 's/CALA/ALA/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CARG/ARG/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CASN/ASN/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CASP/ASP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CCYS/CYS/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CCYX/CYX/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLN/GLN/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLU/GLU/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLY/GLY/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHID/HID/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHIE/HIE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHIP/HIP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHYP/HYP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CILE/ILE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CLEU/LEU/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CLYS/LYS/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CMET/MET/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CPHE/PHE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CPRO/PRO/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CSER/SER/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTHR/THR/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTRP/TRP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTYR/TYR/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CVAL/VAL/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NALA/ALA/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NARG/ARG/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NASN/ASN/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NASP/ASP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NCYS/CYS/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NCYX/CYX/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLN/GLN/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLU/GLU/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLY/GLY/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHID/HID/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHIE/HIE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHIP/HIP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NILE/ILE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NLEU/LEU/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NLYS/LYS/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NMET/MET/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NPHE/PHE/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NPRO/PRO/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NSER/SER/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTHR/THR/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTRP/TRP/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTYR/TYR/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NVAL/VAL/g' ${poi%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/DAS/ASP/g' ${poi%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i 's/CAF/CYS/g' ${poi%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-DIMETHYLARSINOYL-CYSTEINE
sed -i 's/CAS/CYS/g' ${poi%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-(DIMETHYLARSENIC)CYSTEINE
#sed -i 's/AIB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/ALA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  ALA
#sed -i 's/ALM/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/AYA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/BNN/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/CHG/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/CSD/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DAL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DHA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DNP/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/FLA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/HAC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/PRR/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/MAA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/TIH/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/TPQ/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/0CS/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
#sed -i 's/2BU/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    2BU ADE
#sed -i 's/2OP/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    2OP (2S  2-HYDROXYPROPANAL
#sed -i 's/4F3/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    4F3 ALA  CYCLIZED
#sed -i 's/AA4/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
#sed -i 's/ABA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ABA ALA  ALPHA-AMINOBUTYRIC ACID
#sed -i 's/AHO/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
#sed -i 's/AHP/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AHP ALA  2-AMINO-HEPTANOIC ACID
#sed -i 's/AIB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
#sed -i 's/ALA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALA ALA
#sed -i 's/ALC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
#sed -i 's/ALM/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALM ALA  1-METHYL-ALANINAL
#sed -i 's/ALN/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
#sed -i 's/ALS/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
#sed -i 's/ALT/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ALT ALA  THIOALANINE
#sed -i 's/AP7/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AP7 ADE
#sed -i 's/APH/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    APH ALA  P-AMIDINOPHENYL-3-ALANINE
#sed -i 's/AYA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AYA ALA  N-ACETYLALANINE
#sed -i 's/AYG/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    AYG ALA
#sed -i 's/B2A/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    B2A ALA  ALANINE BORONIC ACID
#sed -i 's/B3A/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    B3A ALA  (3S)-3-AMINOBUTANOIC ACID
#sed -i 's/BAL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    BAL ALA  BETA-ALANINE
#sed -i 's/BNN/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
#sed -i 's/C12/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    C12 ALA
#sed -i 's/C99/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    C99 ALA
#sed -i 's/CAB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CAB ALA  4-CARBOXY-4-AMINOBUTANAL
#sed -i 's/CH6/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CH6 ALA
#sed -i 's/CH7/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CH7 ALA
#sed -i 's/CLB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CLB ALA
#sed -i 's/CLD/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CLD ALA
#sed -i 's/CLV/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CLV ALA
#sed -i 's/CQR/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CQR ALA
#sed -i 's/CR2/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CR2 ALA  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/CR5/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CR5 ALA
#sed -i 's/CR7/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CR7 ALA
#sed -i 's/CR8/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CR8 ALA
#sed -i 's/CRK/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CRK ALA
#sed -i 's/CRW/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CRW ALA
#sed -i 's/CRX/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CRX ALA
#sed -i 's/CSI/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CSI ALA
#sed -i 's/CSY/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CSY ALA  MODIFIED TYROSINE COMPLEX
#sed -i 's/CWR/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    CWR ALA
#sed -i 's/DAB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DAB ALA  24-DIAMINOBUTYRIC ACID
#sed -i 's/DAL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DAL ALA  D-ALANINE
#sed -i 's/DAM/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
#sed -i 's/DBU/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
#sed -i 's/DBZ/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
#sed -i 's/DHA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DHA ALA  2-AMINO-ACRYLIC ACID
#sed -i 's/DPP/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    DPP ALA  DIAMMINOPROPANOIC ACID
#sed -i 's/FGL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    FGL ALA  2-AMINOPROPANEDIOIC ACID
#sed -i 's/HHK/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    HHK ALA  (2S)-28-DIAMINOOCTANOIC ACID
#sed -i 's/HMF/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
#sed -i 's/IAM/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
#sed -i 's/IGL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
#sed -i 's/KYN/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    KYN ALA  KYNURENINE
#sed -i 's/LAL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    LAL ALA  NN-DIMETHYL-L-ALANINE
#sed -i 's/MAA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    MAA ALA  N-METHYLALANINE
#sed -i 's/MDO/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    MDO ALA
#sed -i 's/MFC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    MFC ALA  CYCLIZED
#sed -i 's/NAL/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
#sed -i 's/NAM/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    NAM ALA  NAM NAPTHYLAMINOALANINE
#sed -i 's/NCB/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    NCB ALA  CHEMICAL MODIFICATION
#sed -i 's/NRQ/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    NRQ ALA
#sed -i 's/NYC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    NYC ALA
#sed -i 's/ORN/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    ORN ALA  ORNITHINE
#sed -i 's/PIA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    PIA ALA  FUSION OF ALA 65 TYR 66 GLY 67
#sed -i 's/PRR/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
#sed -i 's/PYA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    PYA ALA  3-(110-PHENANTHROL-2-YL)-L-ALANINE
#sed -i 's/PYC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    PYC ALA  PYRROLE-2-CARBOXYLATE
#sed -i 's/PYT/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    PYT ALA  MODIFIED ALANINE
#sed -i 's/RC7/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    RC7 ALA
#sed -i 's/SEC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
#sed -i 's/SIC/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    SIC ALA
#sed -i 's/SUI/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    SUI ALA
#sed -i 's/TIH/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    TIH ALA  BETA(2-THIENYL)ALANINE
#sed -i 's/TPQ/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    TPQ ALA  245-TRIHYDROXYPHENYLALANINE
#sed -i 's/UMA/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    UMA ALA
#sed -i 's/X9Q/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    X9Q ALA
#sed -i 's/XXY/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    XXY ALA
#sed -i 's/XYG/ALA/g' ${poi%.*}_r_prep1.pdb                 ###                    XYG ALA
#sed -i 's/BCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/BUC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/C5C/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/C6C/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CEA/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CME/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSO/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSP/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSX/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSW/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CY1/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CY3/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CYG/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i 's/CYM/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CYS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  CYS
#sed -i 's/CYQ/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/DCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/EFC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/OCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/PEC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/PR3/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCH/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SHC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SMC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SOC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/5CS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    5CS CYS
#sed -i 's/AGT/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    AGT CYS  AGMATINE-CYSTEINE ADDUCT
#sed -i 's/BBC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BBC CYS
#sed -i 's/BCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BCS CYS  BENZYLCYSTEINE
#sed -i 's/BCX/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BCX CYS  BETA-3-CYSTEINE
#sed -i 's/BPE/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BPE CYS
#sed -i 's/BUC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BUC CYS  SS-BUTYLTHIOCYSTEINE
#sed -i 's/C3Y/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    C3Y CYS  MODIFIED CYSTEINE
#sed -i 's/C5C/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
#sed -i 's/C6C/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
#sed -i 's/CAF/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
#sed -i 's/CAS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
#sed -i 's/CCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CCS CYS  CARBOXYMETHYLATED CYSTEINE
#sed -i 's/CME/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CME CYS  MODIFIED CYSTEINE
#sed -i 's/CML/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CML CYS
#sed -i 's/CMT/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CMT CYS  O-METHYLCYSTEINE
#sed -i 's/CS1/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
#sed -i 's/CS3/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CS3 CYS
#sed -i 's/CS4/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CS4 CYS
#sed -i 's/CSA/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSA CYS  S-ACETONYLCYSTEIN
#sed -i 's/CSB/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSB CYS  CYS BOUND TO LEAD ION
#sed -i 's/CSD/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSD CYS  3-SULFINOALANINE
#sed -i 's/CSE/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSE CYS  SELENOCYSTEINE
#sed -i 's/CSO/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSO CYS  INE S-HYDROXYCYSTEINE
#sed -i 's/CSR/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSR CYS  S-ARSONOCYSTEINE
#sed -i 's/CSS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSS CYS  13-THIAZOLE-4-CARBOXYLIC ACID
#sed -i 's/CSU/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSU CYS  CYSTEINE-S-SULFONIC ACID
#sed -i 's/CSW/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSW CYS  CYSTEINE-S-DIOXIDE
#sed -i 's/CSX/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSX CYS  OXOCYSTEINE
#sed -i 's/CSZ/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CSZ CYS  S-SELANYL CYSTEINE
#sed -i 's/CY0/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CY0 CYS  MODIFIED CYSTEINE
#sed -i 's/CY1/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CY1 CYS  ACETAMIDOMETHYLCYSTEINE
#sed -i 's/CY3/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
#sed -i 's/CY4/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CY4 CYS  S-BUTYRYL-CYSTEIN
#sed -i 's/CY7/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CY7 CYS  MODIFIED CYSTEINE
#sed -i 's/CYD/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYD CYS
#sed -i 's/CYF/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
#sed -i 's/CYG/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYG CYS
#sed -i 's/CYQ/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYQ CYS
#sed -i 's/CYR/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYR CYS
#sed -i 's/CYS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYS CYS
#sed -i 's/CZ2/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
#sed -i 's/CZZ/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CZZ CYS  THIARSAHYDROXY-CYSTEINE
#sed -i 's/DCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    DCY CYS  D-CYSTEINE
#sed -i 's/DYS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    DYS CYS
#sed -i 's/EFC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    EFC CYS  SS-(2-FLUOROETHYL)THIOCYSTEINE
#sed -i 's/FOE/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    FOE CYS
#sed -i 's/GT9/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    GT9 CYS  SG ALKYLATED
#sed -i 's/GYC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    GYC CYS
#sed -i 's/HTI/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    HTI CYS
#sed -i 's/KOR/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    KOR CYS  MODIFIED CYSTEINE
#sed -i 's/M0H/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
#sed -i 's/MCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    MCS CYS  MALONYLCYSTEINE
#sed -i 's/NPH/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    NPH CYS
#sed -i 's/NYS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    NYS CYS
#sed -i 's/OCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    OCS CYS  CYSTEINE SULFONIC ACID
#sed -i 's/OCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    OCY CYS  HYDROXYETHYLCYSTEINE
#sed -i 's/P1L/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    P1L CYS  S-PALMITOYL CYSTEINE
#sed -i 's/PBB/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
#sed -i 's/PEC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    PEC CYS  SS-PENTYLTHIOCYSTEINE
#sed -i 's/PR3/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    PR3 CYS  INE DTT-CYSTEINE
#sed -i 's/PYX/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
#sed -i 's/R1A/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    R1A CYS
#sed -i 's/R1B/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    R1B CYS
#sed -i 's/R1F/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    R1F CYS
#sed -i 's/R7A/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    R7A CYS
#sed -i 's/RCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    RCY CYS
#sed -i 's/SAH/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
#sed -i 's/SC2/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SC2 CYS  N-ACETYL-L-CYSTEINE
#sed -i 's/SCH/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SCH CYS  S-METHYL THIOCYSTEINE GROUP
#sed -i 's/SCS/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SCS CYS  MODIFIED CYSTEINE
#sed -i 's/SCY/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SCY CYS  CETYLATED CYSTEINE
#sed -i 's/SHC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SHC CYS  S-HEXYLCYSTEINE
#sed -i 's/SMC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SMC CYS  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/SNC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SNC CYS  S-NITROSO CYSTEINE
#sed -i 's/SOC/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SOC CYS  DIOXYSELENOCYSTEINE
#sed -i 's/TEE/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    TEE CYS  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/TNB/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    TNB CYS  S-(236-TRINITROPHENYL)CYSTEINE
#sed -i 's/TYX/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
#sed -i 's/YCM/CYS/g' ${poi%.*}_r_prep1.pdb                 ###                    YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE
#sed -i 's/2AS/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASA/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASB/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASK/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASL/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASP/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  ASP
#sed -i 's/ASQ/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/BHD/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/DAS/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/DSP/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/3MD/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    3MD ASP  2S3S-3-METHYLASPARTIC ACID
#sed -i 's/A0A/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
#sed -i 's/ACB/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ACB ASP  3-METHYL-ASPARTIC ACID
#sed -i 's/AKL/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
#sed -i 's/ASA/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASA ASP  ASPARTIC ALDEHYDE
#sed -i 's/ASB/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
#sed -i 's/ASI/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASI ASP  L-ISO-ASPARTATE
#sed -i 's/ASK/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
#sed -i 's/ASL/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
#sed -i 's/ASP/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    ASP ASP
#sed -i 's/B3D/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    B3D ASP  3-AMINOPENTANEDIOIC ACID
#sed -i 's/BFD/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
#sed -i 's/BHD/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    BHD ASP  BETA-HYDROXYASPARTIC ACID
#sed -i 's/DAS/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    DAS ASP  D-ASPARTIC ACID
#sed -i 's/DMK/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    DMK ASP  DIMETHYL ASPARTIC ACID
#sed -i 's/IAS/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    IAS ASP  ASPARTYL GROUP
#sed -i 's/OHS/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
#sed -i 's/OXX/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
#sed -i 's/PHD/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
#sed -i 's/SNN/ASP/g' ${poi%.*}_r_prep1.pdb                 ###                    SNN ASP  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/5HP/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/CGU/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/DGL/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/GGL/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/GLU/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  GLU
#sed -i 's/GMA/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/PCA/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/AB7/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    AB7 GLU  ALPHA-AMINOBUTYRIC ACID
#sed -i 's/AR4/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    AR4 GLU
#sed -i 's/B3E/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
#sed -i 's/CGU/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    CGU GLU  CARBOXYLATION OF THE CG ATOM
#sed -i 's/DGL/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    DGL GLU  D-GLU
#sed -i 's/GLU/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    GLU GLU
#sed -i 's/GMA/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    GMA GLU  1-AMIDO-GLUTAMIC ACID
#sed -i 's/ILG/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
#sed -i 's/LME/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
#sed -i 's/MEG/GLU/g' ${poi%.*}_r_prep1.pdb                 ###                    MEG GLU  (2S3R)-3-METHYL-GLUTAMIC ACID
#sed -i 's/DAH/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/DPN/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/HPQ/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/PHE/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  PHE
#sed -i 's/PHI/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/PHL/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/1PA/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
#sed -i 's/23F/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
#sed -i 's/4PH/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    4PH PHE  4-METHYL-L-PHENYLALANINE
#sed -i 's/B2F/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    B2F PHE  PHENYLALANINE BORONIC ACID
#sed -i 's/BIF/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    BIF PHE
#sed -i 's/CHS/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
#sed -i 's/DAH/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    DAH PHE  34-DIHYDROXYDAHNYLALANINE
#sed -i 's/DPH/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    DPH PHE  DEAMINO-METHYL-PHENYLALANINE
#sed -i 's/DPN/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    DPN PHE  D-CONFIGURATION
#sed -i 's/FCL/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    FCL PHE  3-CHLORO-L-PHENYLALANINE
#sed -i 's/FOG/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
#sed -i 's/FRF/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    FRF PHE  PHE FOLLOWED BY REDUCED PHE
#sed -i 's/HPE/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    HPE PHE  HOMOPHENYLALANINE
#sed -i 's/HPH/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    HPH PHE  PHENYLALANINOL GROUP
#sed -i 's/HPQ/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    HPQ PHE  HOMOPHENYLALANINYLMETHANE
#sed -i 's/MEA/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    MEA PHE  N-METHYLPHENYLALANINE
#sed -i 's/MTY/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    MTY PHE  3-HYDROXYPHENYLALANINE
#sed -i 's/NFA/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    NFA PHE  MODIFIED PHENYLALANINE
#sed -i 's/PBF/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
#sed -i 's/PCS/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PCS PHE  PHENYLALANYLMETHYLCHLORIDE
#sed -i 's/PF5/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PF5 PHE  23456-PENTAFLUORO-L-PHENYLALANINE
#sed -i 's/PFF/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PFF PHE  4-FLUORO-L-PHENYLALANINE
#sed -i 's/PHA/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PHA PHE  PHENYLALANINAL
#sed -i 's/PHE/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PHE PHE
#sed -i 's/PHI/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PHI PHE  IODO-PHENYLALANINE
#sed -i 's/PHL/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PHL PHE  L-PHENYLALANINOL
#sed -i 's/PHM/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PHM PHE  PHENYLALANYLMETHANE
#sed -i 's/PM3/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PM3 PHE
#sed -i 's/PPN/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
#sed -i 's/PRQ/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PRQ PHE  PHENYLALANINE
#sed -i 's/PSA/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    PSA PHE
#sed -i 's/SMF/PHE/g' ${poi%.*}_r_prep1.pdb                 ###                    SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
#sed -i 's/GL3/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/GLY/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  GLY
#sed -i 's/GLZ/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/GSC/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/MPQ/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/MSA/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/NMC/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/SAR/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/ACY/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    ACY GLY  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/CHG/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    CHG GLY  CYCLOHEXYL GLYCINE
#sed -i 's/CHP/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
#sed -i 's/GHP/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    GHP GLY  4-HYDROXYPHENYLGLYCINE
#sed -i 's/GL3/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    GL3 GLY  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/GLY/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    GLY GLY
#sed -i 's/GLZ/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    GLZ GLY  AMINO-ACETALDEHYDE
#sed -i 's/GYS/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    GYS GLY
#sed -i 's/IPG/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    IPG GLY  N-ISOPROPYL GLYCINE
#sed -i 's/MEU/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    MEU GLY  O-METHYL-GLYCINE
#sed -i 's/MPQ/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
#sed -i 's/MSA/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    MSA GLY  (2-S-METHYL) SARCOSINE
#sed -i 's/NMC/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
#sed -i 's/PG9/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    PG9 GLY  D-PHENYLGLYCINE
#sed -i 's/SAR/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    SAR GLY  SARCOSINE
#sed -i 's/SHP/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
#sed -i 's/TBG/GLY/g' ${poi%.*}_r_prep1.pdb                 ###                    TBG GLY  T-BUTYL GLYCINE
#sed -i 's/3AH/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/DHI/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HIC/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HIS/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HIS
#sed -i 's/MHS/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/NEM/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/NEP/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HID/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  single delta N protonation
#sed -i 's/HIE/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                  single epsilon N protonation
#sed -i 's/3AH/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    3AH HIS
#sed -i 's/DDE/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    DDE HIS
#sed -i 's/DHI/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    DHI HIS  D-HISTIDINE
#sed -i 's/HIA/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HIA HIS  L-HISTIDINE AMIDE
#sed -i 's/HIC/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HIC HIS  4-METHYL-HISTIDINE
#sed -i 's/HIP/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HIP HIS  ND1-PHOSPHONOHISTIDINE...or commonly used doubly protonated state
#sed -i 's/HIQ/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HIQ HIS  MODIFIED HISTIDINE
#sed -i 's/HIS/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HIS HIS
#sed -i 's/HSO/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    HSO HIS  HISTIDINOL
#sed -i 's/MHS/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    MHS HIS  1-N-METHYLHISTIDINE
#sed -i 's/NEP/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    NEP HIS  N1-PHOSPHONOHISTIDINE
#sed -i 's/NZH/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    NZH HIS
#sed -i 's/OHI/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
#sed -i 's/PSH/HIS/g' ${poi%.*}_r_prep1.pdb                 ###                    PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
#sed -i 's/DIL/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
#sed -i 's/IIL/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
#sed -i 's/ILE/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                  ILE
#sed -i 's/B2I/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    B2I ILE  ISOLEUCINE BORONIC ACID
#sed -i 's/DIL/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    DIL ILE  D-ISOLEUCINE
#sed -i 's/IIL/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    IIL ILE  ISO-ISOLEUCINE
#sed -i 's/ILE/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    ILE ILE
#sed -i 's/ILX/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    ILX ILE  45-DIHYDROXYISOLEUCINE
#sed -i 's/IML/ILE/g' ${poi%.*}_r_prep1.pdb                 ###                    IML ILE  N-METHYLATED
#sed -i 's/ALY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/DLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/KCX/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LLP/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LYM/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LYS/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  LYS
#sed -i 's/LYZ/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/MLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/SHR/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/TRG/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/6CL/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    6CL LYS  6-CARBOXYLYSINE
#sed -i 's/ALY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    ALY LYS  N(6)-ACETYLLYSINE
#sed -i 's/API/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    API LYS  26-DIAMINOPIMELIC ACID
#sed -i 's/APK/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    APK LYS
#sed -i 's/AZK/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
#sed -i 's/B3K/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    B3K LYS  (3S)-37-DIAMINOHEPTANOIC ACID
#sed -i 's/BLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    BLY LYS  LYSINE BORONIC ACID
#sed -i 's/C1X/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    C1X LYS  MODIFIED LYSINE
#sed -i 's/CLG/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CLG LYS
#sed -i 's/CLH/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CLH LYS
#sed -i 's/CYJ/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    CYJ LYS  MODIFIED LYSINE
#sed -i 's/DLS/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    DLS LYS  DI-ACETYL-LYSINE
#sed -i 's/DLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    DLY LYS  D-LYSINE
#sed -i 's/DNL/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    DNL LYS  6-AMINO-HEXANAL
#sed -i 's/FHL/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    FHL LYS  MODIFIED LYSINE
#sed -i 's/GPL/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    GPL LYS  LYSINE GUANOSINE-5-MONOPHOSPHATE
#sed -i 's/IT1/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    IT1 LYS
#sed -i 's/KCX/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    KCX LYS  CARBAMOYLATED LYSINE
#sed -i 's/KGC/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    KGC LYS
#sed -i 's/KST/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
#sed -i 's/LA2/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LA2 LYS
#sed -i 's/LCK/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LCK LYS
#sed -i 's/LCX/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LCX LYS  CARBAMYLATED LYSINE
#sed -i 's/LDH/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LDH LYS  N~6~-ETHYL-L-LYSINE
#sed -i 's/LET/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LET LYS  ODIFIED LYSINE
#sed -i 's/LLP/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LLP LYS
#sed -i 's/LLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
#sed -i 's/LSO/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LSO LYS  MODIFIED LYSINE
#sed -i 's/LYM/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYM LYS  DEOXY-METHYL-LYSINE
#sed -i 's/LYN/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYN LYS  26-DIAMINO-HEXANOIC ACID AMIDE
#sed -i 's/LYP/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
#sed -i 's/LYR/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYR LYS  MODIFIED LYSINE
#sed -i 's/LYS/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYS LYS
#sed -i 's/LYX/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYX LYS  N-(2-COENZYME A)-PROPANOYL-LYSINE
#sed -i 's/LYZ/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    LYZ LYS  5-HYDROXYLYSINE
#sed -i 's/M2L/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    M2L LYS
#sed -i 's/M3L/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    M3L LYS  N-TRIMETHYLLYSINE
#sed -i 's/MCL/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
#sed -i 's/MLY/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    MLY LYS  METHYLATED LYSINE
#sed -i 's/MLZ/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    MLZ LYS  N-METHYL-LYSINE
#sed -i 's/OBS/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    OBS LYS  MODIFIED LYSINE
#sed -i 's/SLZ/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    SLZ LYS  L-THIALYSINE
#sed -i 's/XX1/LYS/g' ${poi%.*}_r_prep1.pdb                 ###                    XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE
#sed -i 's/BUG/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/CLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/DLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/LEU/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  LEU
#sed -i 's/MLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLN/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLP/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/1LU/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
#sed -i 's/2ML/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    2ML LEU  2-METHYLLEUCINE
#sed -i 's/BLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    BLE LEU  LEUCINE BORONIC ACID
#sed -i 's/BUG/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    BUG LEU  TERT-LEUCYL AMINE
#sed -i 's/CLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    CLE LEU  LEUCINE AMIDE
#sed -i 's/DCL/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
#sed -i 's/DLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    DLE LEU  D-LEUCINE
#sed -i 's/DNE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    DNE LEU  D-NORLEUCINE
#sed -i 's/DNG/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    DNG LEU  N-FORMYL-D-NORLEUCINE
#sed -i 's/DNM/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    DNM LEU  D-N-METHYL NORLEUCINE
#sed -i 's/FLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    FLE LEU  FUROYL-LEUCINE
#sed -i 's/HLU/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    HLU LEU  BETA-HYDROXYLEUCINE
#sed -i 's/LED/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    LED LEU  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/LEF/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    LEF LEU  2-5-FLUOROLEUCINE
#sed -i 's/LEU/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    LEU LEU
#sed -i 's/LNT/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    LNT LEU
#sed -i 's/MHL/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    MHL LEU  N-METHYLATED HYDROXY
#sed -i 's/MLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    MLE LEU  N-METHYLATED
#sed -i 's/MLL/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    MLL LEU  METHYL L-LEUCINATE
#sed -i 's/MNL/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    MNL LEU  4N-DIMETHYLNORLEUCINE
#sed -i 's/NLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    NLE LEU  NORLEUCINE
#sed -i 's/NLN/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    NLN LEU  NORLEUCINE AMIDE
#sed -i 's/NLO/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    NLO LEU  O-METHYL-L-NORLEUCINE
#sed -i 's/PLE/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    PLE LEU  LEUCINE PHOSPHINIC ACID
#sed -i 's/PPH/LEU/g' ${poi%.*}_r_prep1.pdb                 ###                    PPH LEU  PHENYLALANINE PHOSPHINIC ACID
#sed -i 's/CXM/MET/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/FME/MET/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/MET/MET/g' ${poi%.*}_r_prep1.pdb                 ###                  MET
#sed -i 's/MSE/MET/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/OMT/MET/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/AME/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    AME MET  ACETYLATED METHIONINE
#sed -i 's/CXM/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    CXM MET  N-CARBOXYMETHIONINE
#sed -i 's/ESC/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
#sed -i 's/FME/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    FME MET  FORMYL-METHIONINE
#sed -i 's/FOR/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    FOR MET
#sed -i 's/MET/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    MET MET
#sed -i 's/MHO/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    MHO MET  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/MME/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    MME MET  N-METHYL METHIONINE
#sed -i 's/MSE/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    MSE MET  ELENOMETHIONINE
#sed -i 's/MSO/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    MSO MET  METHIONINE SULFOXIDE
#sed -i 's/OMT/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    OMT MET  METHIONINE SULFONE
#sed -i 's/SME/MET/g' ${poi%.*}_r_prep1.pdb                 ###                    SME MET  METHIONINE SULFOXIDE
#sed -i 's/ASN/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                  ASN
#sed -i 's/MEN/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASN
#sed -i 's/AFA/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    AFA ASN  N-[7-METHYL-OCT-24-DIENOYL]ASPARAGINE
#sed -i 's/AHB/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    AHB ASN  BETA-HYDROXYASPARAGINE
#sed -i 's/ASN/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    ASN ASN
#sed -i 's/B3X/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    B3X ASN  (3S)-35-DIAMINO-5-OXOPENTANOIC ACID
#sed -i 's/DMH/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    DMH ASN  N4N4-DIMETHYL-ASPARAGINE
#sed -i 's/DSG/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    DSG ASN  D-ASPARAGINE
#sed -i 's/MEN/ASN/g' ${poi%.*}_r_prep1.pdb                 ###                    MEN ASN  GAMMA METHYL ASPARAGINE
#sed -i 's/DPR/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PRO
#sed -i 's/PRO/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                  PRO
#sed -i 's/1AB/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    1AB PRO  14-DIDEOXY-14-IMINO-D-ARABINITOL
#sed -i 's/2MT/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    2MT PRO
#sed -i 's/4FB/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    4FB PRO  (4S)-4-FLUORO-L-PROLINE
#sed -i 's/DPL/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    DPL PRO  4-OXOPROLINE
#sed -i 's/DPR/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    DPR PRO  D-PROLINE
#sed -i 's/H5M/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
#sed -i 's/HY3/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    HY3 PRO  3-HYDROXYPROLINE
#sed -i 's/HYP/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    HYP PRO  4-HYDROXYPROLINE
#sed -i 's/LPD/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    LPD PRO  L-PROLINAMIDE
#sed -i 's/P2Y/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
#sed -i 's/PCA/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    PCA PRO  5-OXOPROLINE
#sed -i 's/POM/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    POM PRO  CIS-5-METHYL-4-OXOPROLINE
#sed -i 's/PRO/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    PRO PRO
#sed -i 's/PRS/PRO/g' ${poi%.*}_r_prep1.pdb                 ###                    PRS PRO  THIOPROLINE
#sed -i 's/DGN/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLN
#sed -i 's/GLN/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                  GLN
#sed -i 's/DGN/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                    DGN GLN  D-GLUTAMINE
#sed -i 's/GHG/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                    GHG GLN  GAMMA-HYDROXY-GLUTAMINE
#sed -i 's/GLH/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                    GLH GLN
#sed -i 's/GLN/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                    GLN GLN
#sed -i 's/MGN/GLN/g' ${poi%.*}_r_prep1.pdb                 ###                    MGN GLN  2-METHYL-GLUTAMINE
#sed -i 's/ACL/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/AGM/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/ARG/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  ARG
#sed -i 's/ARM/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/DAR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/HAR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/HMR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/2MR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    2MR ARG  N3 N4-DIMETHYLARGININE
#sed -i 's/AAR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    AAR ARG  ARGININEAMIDE
#sed -i 's/ACL/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
#sed -i 's/AGM/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    AGM ARG  4-METHYL-ARGININE
#sed -i 's/ALG/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ALG ARG  GUANIDINOBUTYRYL GROUP
#sed -i 's/AR2/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
#sed -i 's/ARG/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ARG ARG
#sed -i 's/ARM/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ARM ARG  DEOXY-METHYL-ARGININE
#sed -i 's/ARO/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ARO ARG  C-GAMMA-HYDROXY ARGININE
#sed -i 's/BOR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    BOR ARG
#sed -i 's/CIR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    CIR ARG  CITRULLINE
#sed -i 's/DA2/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    DA2 ARG  MODIFIED ARGININE
#sed -i 's/DAR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    DAR ARG  D-ARGININE
#sed -i 's/HMR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    HMR ARG  BETA-HOMOARGININE
#sed -i 's/HRG/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    HRG ARG  L-HOMOARGININE
#sed -i 's/MAI/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    MAI ARG  DEOXO-METHYLARGININE
#sed -i 's/MGG/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    MGG ARG  MODIFIED D-ARGININE
#sed -i 's/NMM/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    NMM ARG  MODIFIED ARGININE
#sed -i 's/OPR/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    OPR ARG  C-(3-OXOPROPYL)ARGININE
#sed -i 's/ORQ/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    ORQ ARG  N~5~-ACETYL-L-ORNITHINE
#sed -i 's/TYZ/ARG/g' ${poi%.*}_r_prep1.pdb                 ###                    TYZ ARG  PARA ACETAMIDO BENZOIC ACID
#sed -i 's/DSN/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/MIS/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/OAS/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SAC/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SEL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SEP/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SER/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  SER
#sed -i 's/SET/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SVA/SER/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/B3S/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
#sed -i 's/BG1/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    BG1 SER
#sed -i 's/DHL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    DHL SER  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/DSE/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    DSE SER  D-SERINE N-METHYLATED
#sed -i 's/DSN/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    DSN SER  D-SERINE
#sed -i 's/FGP/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    FGP SER
#sed -i 's/GVL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
#sed -i 's/HSE/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    HSE SER  L-HOMOSERINE
#sed -i 's/HSL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    HSL SER  HOMOSERINE LACTONE
#sed -i 's/MC1/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    MC1 SER  METHICILLIN ACYL-SERINE
#sed -i 's/MIS/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    MIS SER  MODIFIED SERINE
#sed -i 's/N10/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
#sed -i 's/NC1/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    NC1 SER  NITROCEFIN ACYL-SERINE
#sed -i 's/OAS/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    OAS SER  O-ACETYLSERINE
#sed -i 's/OSE/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    OSE SER  O-SULFO-L-SERINE
#sed -i 's/PG1/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
#sed -i 's/PYR/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    PYR SER  CHEMICALLY MODIFIED
#sed -i 's/S1H/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
#sed -i 's/SAC/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SAC SER  N-ACETYL-SERINE
#sed -i 's/SBD/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SBD SER
#sed -i 's/SBG/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SBG SER  MODIFIED SERINE
#sed -i 's/SBL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SBL SER
#sed -i 's/SDP/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SDP SER
#sed -i 's/SEB/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SEB SER  O-BENZYLSULFONYL-SERINE
#sed -i 's/SEL/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SEL SER  2-AMINO-13-PROPANEDIOL
#sed -i 's/SEP/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SEP SER  E PHOSPHOSERINE
#sed -i 's/SER/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SER SER
#sed -i 's/SET/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SET SER  AMINOSERINE
#sed -i 's/SGB/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SGB SER  MODIFIED SERINE
#sed -i 's/SGR/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SGR SER  MODIFIED SERINE
#sed -i 's/SOY/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SOY SER  OXACILLOYL-ACYLATED SERINE
#sed -i 's/SUN/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SUN SER  TABUN CONJUGATED SERINE
#sed -i 's/SVA/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SVA SER  SERINE VANADATE
#sed -i 's/SVV/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SVV SER  MODIFIED SERINE
#sed -i 's/SVX/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SVX SER  MODIFIED SERINE
#sed -i 's/SVY/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SVY SER  MODIFIED SERINE
#sed -i 's/SVZ/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SVZ SER  MODIFIED SERINE
#sed -i 's/SXE/SER/g' ${poi%.*}_r_prep1.pdb                 ###                    SXE SER  MODIFIED SERINE
#sed -i 's/ALO/THR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/BMT/THR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/DTH/THR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/THR/THR/g' ${poi%.*}_r_prep1.pdb                 ###                  THR
#sed -i 's/TPO/THR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/AEI/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    AEI THR  ACYLATED THR
#sed -i 's/ALO/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    ALO THR  ALLO-THREONINE
#sed -i 's/BMT/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    BMT THR
#sed -i 's/CRO/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    CRO THR  CYCLIZED
#sed -i 's/CTH/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    CTH THR  4-CHLOROTHREONINE
#sed -i 's/DTH/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    DTH THR  D-THREONINE
#sed -i 's/OLT/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    OLT THR  O-METHYL-L-THREONINE
#sed -i 's/TBM/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    TBM THR
#sed -i 's/TH5/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    TH5 THR  O-ACETYL-L-THREONINE
#sed -i 's/THC/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    THC THR  N-METHYLCARBONYLTHREONINE
#sed -i 's/THR/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    THR THR
#sed -i 's/TMD/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    TMD THR  N-METHYLATED EPSILON C ALKYLATED
#sed -i 's/TPO/THR/g' ${poi%.*}_r_prep1.pdb                 ###                    TPO THR  HOSPHOTHREONINE
#sed -i 's/DIV/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/DVA/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/MVA/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/VAL/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                  VAL
#sed -i 's/B2V/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    B2V VAL  VALINE BORONIC ACID
#sed -i 's/DIV/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    DIV VAL  D-ISOVALINE
#sed -i 's/DVA/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    DVA VAL  D-VALINE
#sed -i 's/MNV/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    MNV VAL  N-METHYL-C-AMINO VALINE
#sed -i 's/MVA/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    MVA VAL  N-METHYLATED
#sed -i 's/NVA/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    NVA VAL  NORVALINE
#sed -i 's/VAD/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VAD VAL  DEAMINOHYDROXYVALINE
#sed -i 's/VAF/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VAF VAL  METHYLVALINE
#sed -i 's/VAL/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VAL VAL
#sed -i 's/VDL/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VDL VAL  (2R3R)-23-DIAMINOBUTANOIC ACID
#sed -i 's/VLL/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VLL VAL  (2S)-23-DIAMINOBUTANOIC ACID
#sed -i 's/VME/VAL/g' ${poi%.*}_r_prep1.pdb                 ###                    VME VAL  O- METHYLVALINE
#sed -i 's/DTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/HTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/LTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TPL/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TRO/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TRP/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                  TRP
#sed -i 's/BTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    BTR TRP  6-BROMO-TRYPTOPHAN
#sed -i 's/1TQ/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
#sed -i 's/23S/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    23S TRP  MODIFIED TRYPTOPHAN
#sed -i 's/32S/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    32S TRP  MODIFIED TRYPTOPHAN
#sed -i 's/32T/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    32T TRP  MODIFIED TRYPTOPHAN
#sed -i 's/4DP/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    4DP TRP
#sed -i 's/4FW/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    4FW TRP  4-FLUOROTRYPTOPHANE
#sed -i 's/4HT/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    4HT TRP  4-HYDROXYTRYPTOPHAN
#sed -i 's/4IN/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    4IN TRP  4-AMINO-L-TRYPTOPHAN
#sed -i 's/6CW/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    6CW TRP  6-CHLORO-L-TRYPTOPHAN
#sed -i 's/DTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    DTR TRP  D-TRYPTOPHAN
#sed -i 's/FTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    FTR TRP  FLUOROTRYPTOPHANE
#sed -i 's/HTR/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    HTR TRP  BETA-HYDROXYTRYPTOPHANE
#sed -i 's/PAT/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
#sed -i 's/TOX/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TOX TRP
#sed -i 's/TPL/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TPL TRP  TRYTOPHANOL
#sed -i 's/TQQ/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TQQ TRP
#sed -i 's/TRF/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRF TRP  N1-FORMYL-TRYPTOPHAN
#sed -i 's/TRN/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRN TRP  AZA-TRYPTOPHAN
#sed -i 's/TRO/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRO TRP  2-HYDROXY-TRYPTOPHAN
#sed -i 's/TRP/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRP TRP
#sed -i 's/TRQ/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRQ TRP
#sed -i 's/TRW/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRW TRP
#sed -i 's/TRX/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TRX TRP  6-HYDROXYTRYPTOPHAN
#sed -i 's/TTQ/TRP/g' ${poi%.*}_r_prep1.pdb                 ###                    TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
#sed -i 's/DTY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/IYR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/PAQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/PTR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/STY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYB/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYS/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  TYR
#sed -i 's/TYY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/1TY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    1TY TYR
#sed -i 's/2TY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    2TY TYR
#sed -i 's/3TY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    3TY TYR  MODIFIED TYROSINE
#sed -i 's/B3Y/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    B3Y TYR
#sed -i 's/CRQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    CRQ TYR
#sed -i 's/DBY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    DBY TYR  35 DIBROMOTYROSINE
#sed -i 's/DPQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    DPQ TYR  TYROSINE DERIVATIVE
#sed -i 's/DTY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    DTY TYR  D-TYROSINE
#sed -i 's/ESB/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    ESB TYR
#sed -i 's/FLT/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    FLT TYR  FLUOROMALONYL TYROSINE
#sed -i 's/FTY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
#sed -i 's/IYR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    IYR TYR  3-IODO-TYROSINE
#sed -i 's/MBQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    MBQ TYR
#sed -i 's/NIY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    NIY TYR  META-NITRO-TYROSINE
#sed -i 's/NBQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    NBQ TYR
#sed -i 's/OTY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    OTY TYR
#sed -i 's/PAQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    PAQ TYR  SEE REMARK 999
#sed -i 's/PTH/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
#sed -i 's/PTM/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
#sed -i 's/PTR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    PTR TYR  O-PHOSPHOTYROSINE
#sed -i 's/TCQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TCQ TYR  MODIFIED TYROSINE
#sed -i 's/TTS/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TTS TYR
#sed -i 's/TY2/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TY2 TYR  3-AMINO-L-TYROSINE
#sed -i 's/TY3/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TY3 TYR  3-HYDROXY-L-TYROSINE
#sed -i 's/TYB/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYB TYR  TYROSINAL
#sed -i 's/TYC/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYC TYR  L-TYROSINAMIDE
#sed -i 's/TYI/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYI TYR  35-DIIODOTYROSINE
#sed -i 's/TYN/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYN TYR  ADDUCT AT HYDROXY GROUP
#sed -i 's/TYO/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYO TYR
#sed -i 's/TYQ/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
#sed -i 's/TYR/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYR TYR
#sed -i 's/TYS/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYS TYR  INE SULPHONATED TYROSINE
#sed -i 's/TYT/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYT TYR
#sed -i 's/TYY/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
#sed -i 's/YOF/TYR/g' ${poi%.*}_r_prep1.pdb                 ###                    YOF TYR  3-FLUOROTYROSINE


#################################################################################################################################################################

sed -i 's/CSD/CYS/g' ${e3%.*}_r_prep1.pdb                 #3-SULFINOALANINE
sed -i 's/HYP/PRO/g' ${e3%.*}_r_prep1.pdb                 #4-HYDROXYPROLINE
sed -i 's/BMT/THR/g' ${e3%.*}_r_prep1.pdb                 #4-METHYL-4-[(E)-2-BUTENYL]-4,N-METHYL-THREONINE
sed -i 's/5HP/GLU/g' ${e3%.*}_r_prep1.pdb                 #5-HYDROXYPROLINE
sed -i 's/ABA/ALA/g' ${e3%.*}_r_prep1.pdb                 #ALPHA-AMINOBUTYRIC_ACID
sed -i 's/AIB/ALA/g' ${e3%.*}_r_prep1.pdb                 #ALPHA-AMINOISOBUTYRIC_ACID
sed -i 's/CSW/CYS/g' ${e3%.*}_r_prep1.pdb                 #CYSTEINE-S-DIOXIDE
sed -i 's/OCS/CYS/g' ${e3%.*}_r_prep1.pdb                 #CYSTEINESULFONIC_ACID
sed -i 's/DAL/ALA/g' ${e3%.*}_r_prep1.pdb                 #D-ALANINE
sed -i 's/DAR/ARG/g' ${e3%.*}_r_prep1.pdb                 #D-ARGININE
sed -i 's/DSG/ASN/g' ${e3%.*}_r_prep1.pdb                 #D-ASPARAGINE
sed -i 's/DSP/ASP/g' ${e3%.*}_r_prep1.pdb                 #D-ASPARTATE
sed -i 's/DCY/CYS/g' ${e3%.*}_r_prep1.pdb                 #D-CYSTEINE
sed -i 's/CRO/CRO/g' ${e3%.*}_r_prep1.pdb                 #DECARBOXY(PARAHYDROXYBENZYLIDENE-IMIDAZOLIDINONE)THREONINE
sed -i 's/DGL/GLU/g' ${e3%.*}_r_prep1.pdb                 #D-GLUTAMATE
sed -i 's/DGN/GLN/g' ${e3%.*}_r_prep1.pdb                 #D-GLUTAMINE
sed -i 's/DHI/HIS/g' ${e3%.*}_r_prep1.pdb                 #D-HISTIDINE
sed -i 's/DIL/ILE/g' ${e3%.*}_r_prep1.pdb                 #D-ISOLEUCINE
sed -i 's/DIV/VAL/g' ${e3%.*}_r_prep1.pdb                 #D-ISOVALINE
sed -i 's/DLE/LEU/g' ${e3%.*}_r_prep1.pdb                 #D-LEUCINE
sed -i 's/DLY/LYS/g' ${e3%.*}_r_prep1.pdb                 #D-LYSINE
sed -i 's/DPN/PHE/g' ${e3%.*}_r_prep1.pdb                 #D-PHENYLALANINE
sed -i 's/DPR/PRO/g' ${e3%.*}_r_prep1.pdb                 #D-PROLINE
sed -i 's/DSN/SER/g' ${e3%.*}_r_prep1.pdb                 #D-SERINE
sed -i 's/DTH/THR/g' ${e3%.*}_r_prep1.pdb                 #D-THREONINE
sed -i 's/DTR/DTR/g' ${e3%.*}_r_prep1.pdb                 #D-TRYPTOPHANE
sed -i 's/DTY/TYR/g' ${e3%.*}_r_prep1.pdb                 #D-TYROSINE
sed -i 's/DVA/VAL/g' ${e3%.*}_r_prep1.pdb                 #D-VALINE
sed -i 's/CGU/GLU/g' ${e3%.*}_r_prep1.pdb                 #GAMMA-CARBOXY-GLUTAMIC_ACID
sed -i 's/KCX/LYS/g' ${e3%.*}_r_prep1.pdb                 #LYSINE_NZ-CARBOXYLIC_ACID
sed -i 's/LLP/LYS/g' ${e3%.*}_r_prep1.pdb                 #LYSINE-PYRIDOXAL-5'-PHOSPHATE
sed -i 's/CXM/MET/g' ${e3%.*}_r_prep1.pdb                 #N-CARBOXYMETHIONINE
sed -i 's/FME/MET/g' ${e3%.*}_r_prep1.pdb                 #N-FORMYLMETHIONINE
sed -i 's/MLE/LEU/g' ${e3%.*}_r_prep1.pdb                 #N-METHYLLEUCINE
sed -i 's/MVA/VAL/g' ${e3%.*}_r_prep1.pdb                 #N-METHYLVALINE
sed -i 's/NLE/LEU/g' ${e3%.*}_r_prep1.pdb                 #NORLEUCINE
sed -i 's/PTR/TYR/g' ${e3%.*}_r_prep1.pdb                 #O-PHOSPHOTYROSINE
sed -i 's/ORN/ALA/g' ${e3%.*}_r_prep1.pdb                 #ORNITHINE
sed -i 's/SEP/SER/g' ${e3%.*}_r_prep1.pdb                 #PHOSPHOSERINE
sed -i 's/TPO/THR/g' ${e3%.*}_r_prep1.pdb                 #PHOSPHOTHREONINE
sed -i 's/PCA/GLU/g' ${e3%.*}_r_prep1.pdb                 #PYROGLUTAMIC_ACID
sed -i 's/SAR/GLY/g' ${e3%.*}_r_prep1.pdb                 #SARCOSINE
sed -i 's/CEA/CYS/g' ${e3%.*}_r_prep1.pdb                 #S-HYDROXY-CYSTEINE
sed -i 's/CSO/CYS/g' ${e3%.*}_r_prep1.pdb                 #S-HYDROXYCYSTEINE
sed -i 's/CSS/CYS/g' ${e3%.*}_r_prep1.pdb                 #S-MERCAPTOCYSTEINE
sed -i 's/CSX/CYS/g' ${e3%.*}_r_prep1.pdb                 #S-OXY_CYSTEINE
sed -i 's/CME/CYS/g' ${e3%.*}_r_prep1.pdb                 #S,S-(2-HYDROXYETHYL)THIOCYSTEINE
sed -i 's/TYS/TYR/g' ${e3%.*}_r_prep1.pdb                 #SULFONATED_TYROSINE
sed -i 's/TPQ/PHE/g' ${e3%.*}_r_prep1.pdb                 #TOPO-QUINONE
sed -i 's/STY/TYR/g' ${e3%.*}_r_prep1.pdb                 #TYROSINE-O-SULPHONIC_ACID
sed -i 's/CCS/CYS/g' ${e3%.*}_r_prep1.pdb                 #https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i 's/CALA/ALA/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CARG/ARG/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CASN/ASN/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CASP/ASP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CCYS/CYS/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CCYX/CYX/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLN/GLN/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLU/GLU/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLY/GLY/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHID/HID/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHIE/HIE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHIP/HIP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHYP/HYP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CILE/ILE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CLEU/LEU/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CLYS/LYS/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CMET/MET/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CPHE/PHE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CPRO/PRO/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CSER/SER/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTHR/THR/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTRP/TRP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTYR/TYR/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CVAL/VAL/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NALA/ALA/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NARG/ARG/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NASN/ASN/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NASP/ASP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NCYS/CYS/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NCYX/CYX/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLN/GLN/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLU/GLU/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLY/GLY/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHID/HID/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHIE/HIE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHIP/HIP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NILE/ILE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NLEU/LEU/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NLYS/LYS/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NMET/MET/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NPHE/PHE/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NPRO/PRO/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NSER/SER/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTHR/THR/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTRP/TRP/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTYR/TYR/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NVAL/VAL/g' ${e3%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/DAS/ASP/g' ${e3%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i 's/CAF/CYS/g' ${e3%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-DIMETHYLARSINOYL-CYSTEINE
sed -i 's/CAS/CYS/g' ${e3%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-(DIMETHYLARSENIC)CYSTEINE
#sed -i 's/AIB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/ALA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  ALA
#sed -i 's/ALM/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/AYA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/BNN/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/CHG/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/CSD/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DAL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DHA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DNP/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/FLA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/HAC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/PRR/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/MAA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/TIH/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/TPQ/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/0CS/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
#sed -i 's/2BU/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    2BU ADE
#sed -i 's/2OP/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    2OP (2S  2-HYDROXYPROPANAL
#sed -i 's/4F3/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    4F3 ALA  CYCLIZED
#sed -i 's/AA4/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
#sed -i 's/ABA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ABA ALA  ALPHA-AMINOBUTYRIC ACID
#sed -i 's/AHO/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
#sed -i 's/AHP/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AHP ALA  2-AMINO-HEPTANOIC ACID
#sed -i 's/AIB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
#sed -i 's/ALA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALA ALA
#sed -i 's/ALC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
#sed -i 's/ALM/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALM ALA  1-METHYL-ALANINAL
#sed -i 's/ALN/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
#sed -i 's/ALS/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
#sed -i 's/ALT/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ALT ALA  THIOALANINE
#sed -i 's/AP7/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AP7 ADE
#sed -i 's/APH/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    APH ALA  P-AMIDINOPHENYL-3-ALANINE
#sed -i 's/AYA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AYA ALA  N-ACETYLALANINE
#sed -i 's/AYG/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    AYG ALA
#sed -i 's/B2A/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    B2A ALA  ALANINE BORONIC ACID
#sed -i 's/B3A/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    B3A ALA  (3S)-3-AMINOBUTANOIC ACID
#sed -i 's/BAL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    BAL ALA  BETA-ALANINE
#sed -i 's/BNN/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
#sed -i 's/C12/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    C12 ALA
#sed -i 's/C99/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    C99 ALA
#sed -i 's/CAB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CAB ALA  4-CARBOXY-4-AMINOBUTANAL
#sed -i 's/CH6/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CH6 ALA
#sed -i 's/CH7/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CH7 ALA
#sed -i 's/CLB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CLB ALA
#sed -i 's/CLD/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CLD ALA
#sed -i 's/CLV/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CLV ALA
#sed -i 's/CQR/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CQR ALA
#sed -i 's/CR2/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CR2 ALA  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/CR5/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CR5 ALA
#sed -i 's/CR7/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CR7 ALA
#sed -i 's/CR8/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CR8 ALA
#sed -i 's/CRK/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CRK ALA
#sed -i 's/CRW/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CRW ALA
#sed -i 's/CRX/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CRX ALA
#sed -i 's/CSI/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CSI ALA
#sed -i 's/CSY/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CSY ALA  MODIFIED TYROSINE COMPLEX
#sed -i 's/CWR/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    CWR ALA
#sed -i 's/DAB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DAB ALA  24-DIAMINOBUTYRIC ACID
#sed -i 's/DAL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DAL ALA  D-ALANINE
#sed -i 's/DAM/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
#sed -i 's/DBU/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
#sed -i 's/DBZ/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
#sed -i 's/DHA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DHA ALA  2-AMINO-ACRYLIC ACID
#sed -i 's/DPP/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    DPP ALA  DIAMMINOPROPANOIC ACID
#sed -i 's/FGL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    FGL ALA  2-AMINOPROPANEDIOIC ACID
#sed -i 's/HHK/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    HHK ALA  (2S)-28-DIAMINOOCTANOIC ACID
#sed -i 's/HMF/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
#sed -i 's/IAM/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
#sed -i 's/IGL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
#sed -i 's/KYN/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    KYN ALA  KYNURENINE
#sed -i 's/LAL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    LAL ALA  NN-DIMETHYL-L-ALANINE
#sed -i 's/MAA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    MAA ALA  N-METHYLALANINE
#sed -i 's/MDO/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    MDO ALA
#sed -i 's/MFC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    MFC ALA  CYCLIZED
#sed -i 's/NAL/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
#sed -i 's/NAM/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    NAM ALA  NAM NAPTHYLAMINOALANINE
#sed -i 's/NCB/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    NCB ALA  CHEMICAL MODIFICATION
#sed -i 's/NRQ/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    NRQ ALA
#sed -i 's/NYC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    NYC ALA
#sed -i 's/ORN/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    ORN ALA  ORNITHINE
#sed -i 's/PIA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    PIA ALA  FUSION OF ALA 65 TYR 66 GLY 67
#sed -i 's/PRR/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
#sed -i 's/PYA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    PYA ALA  3-(110-PHENANTHROL-2-YL)-L-ALANINE
#sed -i 's/PYC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    PYC ALA  PYRROLE-2-CARBOXYLATE
#sed -i 's/PYT/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    PYT ALA  MODIFIED ALANINE
#sed -i 's/RC7/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    RC7 ALA
#sed -i 's/SEC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
#sed -i 's/SIC/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    SIC ALA
#sed -i 's/SUI/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    SUI ALA
#sed -i 's/TIH/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    TIH ALA  BETA(2-THIENYL)ALANINE
#sed -i 's/TPQ/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    TPQ ALA  245-TRIHYDROXYPHENYLALANINE
#sed -i 's/UMA/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    UMA ALA
#sed -i 's/X9Q/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    X9Q ALA
#sed -i 's/XXY/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    XXY ALA
#sed -i 's/XYG/ALA/g' ${e3%.*}_r_prep1.pdb                 ###                    XYG ALA
#sed -i 's/BCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/BUC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/C5C/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/C6C/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CEA/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CME/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSO/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSP/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSX/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSW/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CY1/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CY3/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CYG/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i 's/CYM/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CYS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  CYS
#sed -i 's/CYQ/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/DCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/EFC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/OCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/PEC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/PR3/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCH/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SHC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SMC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SOC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/5CS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    5CS CYS
#sed -i 's/AGT/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    AGT CYS  AGMATINE-CYSTEINE ADDUCT
#sed -i 's/BBC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BBC CYS
#sed -i 's/BCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BCS CYS  BENZYLCYSTEINE
#sed -i 's/BCX/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BCX CYS  BETA-3-CYSTEINE
#sed -i 's/BPE/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BPE CYS
#sed -i 's/BUC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BUC CYS  SS-BUTYLTHIOCYSTEINE
#sed -i 's/C3Y/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    C3Y CYS  MODIFIED CYSTEINE
#sed -i 's/C5C/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
#sed -i 's/C6C/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
#sed -i 's/CAF/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
#sed -i 's/CAS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
#sed -i 's/CCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CCS CYS  CARBOXYMETHYLATED CYSTEINE
#sed -i 's/CME/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CME CYS  MODIFIED CYSTEINE
#sed -i 's/CML/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CML CYS
#sed -i 's/CMT/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CMT CYS  O-METHYLCYSTEINE
#sed -i 's/CS1/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
#sed -i 's/CS3/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CS3 CYS
#sed -i 's/CS4/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CS4 CYS
#sed -i 's/CSA/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSA CYS  S-ACETONYLCYSTEIN
#sed -i 's/CSB/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSB CYS  CYS BOUND TO LEAD ION
#sed -i 's/CSD/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSD CYS  3-SULFINOALANINE
#sed -i 's/CSE/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSE CYS  SELENOCYSTEINE
#sed -i 's/CSO/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSO CYS  INE S-HYDROXYCYSTEINE
#sed -i 's/CSR/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSR CYS  S-ARSONOCYSTEINE
#sed -i 's/CSS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSS CYS  13-THIAZOLE-4-CARBOXYLIC ACID
#sed -i 's/CSU/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSU CYS  CYSTEINE-S-SULFONIC ACID
#sed -i 's/CSW/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSW CYS  CYSTEINE-S-DIOXIDE
#sed -i 's/CSX/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSX CYS  OXOCYSTEINE
#sed -i 's/CSZ/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CSZ CYS  S-SELANYL CYSTEINE
#sed -i 's/CY0/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CY0 CYS  MODIFIED CYSTEINE
#sed -i 's/CY1/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CY1 CYS  ACETAMIDOMETHYLCYSTEINE
#sed -i 's/CY3/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
#sed -i 's/CY4/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CY4 CYS  S-BUTYRYL-CYSTEIN
#sed -i 's/CY7/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CY7 CYS  MODIFIED CYSTEINE
#sed -i 's/CYD/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYD CYS
#sed -i 's/CYF/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
#sed -i 's/CYG/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYG CYS
#sed -i 's/CYQ/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYQ CYS
#sed -i 's/CYR/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYR CYS
#sed -i 's/CYS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYS CYS
#sed -i 's/CZ2/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
#sed -i 's/CZZ/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CZZ CYS  THIARSAHYDROXY-CYSTEINE
#sed -i 's/DCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    DCY CYS  D-CYSTEINE
#sed -i 's/DYS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    DYS CYS
#sed -i 's/EFC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    EFC CYS  SS-(2-FLUOROETHYL)THIOCYSTEINE
#sed -i 's/FOE/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    FOE CYS
#sed -i 's/GT9/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    GT9 CYS  SG ALKYLATED
#sed -i 's/GYC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    GYC CYS
#sed -i 's/HTI/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    HTI CYS
#sed -i 's/KOR/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    KOR CYS  MODIFIED CYSTEINE
#sed -i 's/M0H/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
#sed -i 's/MCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    MCS CYS  MALONYLCYSTEINE
#sed -i 's/NPH/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    NPH CYS
#sed -i 's/NYS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    NYS CYS
#sed -i 's/OCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    OCS CYS  CYSTEINE SULFONIC ACID
#sed -i 's/OCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    OCY CYS  HYDROXYETHYLCYSTEINE
#sed -i 's/P1L/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    P1L CYS  S-PALMITOYL CYSTEINE
#sed -i 's/PBB/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
#sed -i 's/PEC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    PEC CYS  SS-PENTYLTHIOCYSTEINE
#sed -i 's/PR3/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    PR3 CYS  INE DTT-CYSTEINE
#sed -i 's/PYX/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
#sed -i 's/R1A/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    R1A CYS
#sed -i 's/R1B/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    R1B CYS
#sed -i 's/R1F/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    R1F CYS
#sed -i 's/R7A/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    R7A CYS
#sed -i 's/RCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    RCY CYS
#sed -i 's/SAH/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
#sed -i 's/SC2/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SC2 CYS  N-ACETYL-L-CYSTEINE
#sed -i 's/SCH/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SCH CYS  S-METHYL THIOCYSTEINE GROUP
#sed -i 's/SCS/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SCS CYS  MODIFIED CYSTEINE
#sed -i 's/SCY/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SCY CYS  CETYLATED CYSTEINE
#sed -i 's/SHC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SHC CYS  S-HEXYLCYSTEINE
#sed -i 's/SMC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SMC CYS  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/SNC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SNC CYS  S-NITROSO CYSTEINE
#sed -i 's/SOC/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SOC CYS  DIOXYSELENOCYSTEINE
#sed -i 's/TEE/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    TEE CYS  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/TNB/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    TNB CYS  S-(236-TRINITROPHENYL)CYSTEINE
#sed -i 's/TYX/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
#sed -i 's/YCM/CYS/g' ${e3%.*}_r_prep1.pdb                 ###                    YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE
#sed -i 's/2AS/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASA/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASB/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASK/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASL/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASP/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  ASP
#sed -i 's/ASQ/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/BHD/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/DAS/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/DSP/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/3MD/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    3MD ASP  2S3S-3-METHYLASPARTIC ACID
#sed -i 's/A0A/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
#sed -i 's/ACB/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ACB ASP  3-METHYL-ASPARTIC ACID
#sed -i 's/AKL/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
#sed -i 's/ASA/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASA ASP  ASPARTIC ALDEHYDE
#sed -i 's/ASB/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
#sed -i 's/ASI/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASI ASP  L-ISO-ASPARTATE
#sed -i 's/ASK/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
#sed -i 's/ASL/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
#sed -i 's/ASP/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    ASP ASP
#sed -i 's/B3D/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    B3D ASP  3-AMINOPENTANEDIOIC ACID
#sed -i 's/BFD/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
#sed -i 's/BHD/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    BHD ASP  BETA-HYDROXYASPARTIC ACID
#sed -i 's/DAS/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    DAS ASP  D-ASPARTIC ACID
#sed -i 's/DMK/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    DMK ASP  DIMETHYL ASPARTIC ACID
#sed -i 's/IAS/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    IAS ASP  ASPARTYL GROUP
#sed -i 's/OHS/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
#sed -i 's/OXX/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
#sed -i 's/PHD/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
#sed -i 's/SNN/ASP/g' ${e3%.*}_r_prep1.pdb                 ###                    SNN ASP  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/5HP/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/CGU/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/DGL/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/GGL/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/GLU/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  GLU
#sed -i 's/GMA/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/PCA/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/AB7/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    AB7 GLU  ALPHA-AMINOBUTYRIC ACID
#sed -i 's/AR4/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    AR4 GLU
#sed -i 's/B3E/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
#sed -i 's/CGU/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    CGU GLU  CARBOXYLATION OF THE CG ATOM
#sed -i 's/DGL/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    DGL GLU  D-GLU
#sed -i 's/GLU/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    GLU GLU
#sed -i 's/GMA/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    GMA GLU  1-AMIDO-GLUTAMIC ACID
#sed -i 's/ILG/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
#sed -i 's/LME/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
#sed -i 's/MEG/GLU/g' ${e3%.*}_r_prep1.pdb                 ###                    MEG GLU  (2S3R)-3-METHYL-GLUTAMIC ACID
#sed -i 's/DAH/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/DPN/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/HPQ/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/PHE/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  PHE
#sed -i 's/PHI/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/PHL/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/1PA/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
#sed -i 's/23F/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
#sed -i 's/4PH/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    4PH PHE  4-METHYL-L-PHENYLALANINE
#sed -i 's/B2F/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    B2F PHE  PHENYLALANINE BORONIC ACID
#sed -i 's/BIF/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    BIF PHE
#sed -i 's/CHS/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
#sed -i 's/DAH/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    DAH PHE  34-DIHYDROXYDAHNYLALANINE
#sed -i 's/DPH/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    DPH PHE  DEAMINO-METHYL-PHENYLALANINE
#sed -i 's/DPN/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    DPN PHE  D-CONFIGURATION
#sed -i 's/FCL/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    FCL PHE  3-CHLORO-L-PHENYLALANINE
#sed -i 's/FOG/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
#sed -i 's/FRF/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    FRF PHE  PHE FOLLOWED BY REDUCED PHE
#sed -i 's/HPE/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    HPE PHE  HOMOPHENYLALANINE
#sed -i 's/HPH/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    HPH PHE  PHENYLALANINOL GROUP
#sed -i 's/HPQ/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    HPQ PHE  HOMOPHENYLALANINYLMETHANE
#sed -i 's/MEA/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    MEA PHE  N-METHYLPHENYLALANINE
#sed -i 's/MTY/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    MTY PHE  3-HYDROXYPHENYLALANINE
#sed -i 's/NFA/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    NFA PHE  MODIFIED PHENYLALANINE
#sed -i 's/PBF/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
#sed -i 's/PCS/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PCS PHE  PHENYLALANYLMETHYLCHLORIDE
#sed -i 's/PF5/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PF5 PHE  23456-PENTAFLUORO-L-PHENYLALANINE
#sed -i 's/PFF/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PFF PHE  4-FLUORO-L-PHENYLALANINE
#sed -i 's/PHA/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PHA PHE  PHENYLALANINAL
#sed -i 's/PHE/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PHE PHE
#sed -i 's/PHI/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PHI PHE  IODO-PHENYLALANINE
#sed -i 's/PHL/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PHL PHE  L-PHENYLALANINOL
#sed -i 's/PHM/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PHM PHE  PHENYLALANYLMETHANE
#sed -i 's/PM3/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PM3 PHE
#sed -i 's/PPN/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
#sed -i 's/PRQ/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PRQ PHE  PHENYLALANINE
#sed -i 's/PSA/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    PSA PHE
#sed -i 's/SMF/PHE/g' ${e3%.*}_r_prep1.pdb                 ###                    SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
#sed -i 's/GL3/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/GLY/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  GLY
#sed -i 's/GLZ/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/GSC/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/MPQ/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/MSA/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/NMC/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/SAR/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/ACY/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    ACY GLY  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/CHG/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    CHG GLY  CYCLOHEXYL GLYCINE
#sed -i 's/CHP/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
#sed -i 's/GHP/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    GHP GLY  4-HYDROXYPHENYLGLYCINE
#sed -i 's/GL3/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    GL3 GLY  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/GLY/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    GLY GLY
#sed -i 's/GLZ/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    GLZ GLY  AMINO-ACETALDEHYDE
#sed -i 's/GYS/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    GYS GLY
#sed -i 's/IPG/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    IPG GLY  N-ISOPROPYL GLYCINE
#sed -i 's/MEU/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    MEU GLY  O-METHYL-GLYCINE
#sed -i 's/MPQ/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
#sed -i 's/MSA/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    MSA GLY  (2-S-METHYL) SARCOSINE
#sed -i 's/NMC/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
#sed -i 's/PG9/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    PG9 GLY  D-PHENYLGLYCINE
#sed -i 's/SAR/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    SAR GLY  SARCOSINE
#sed -i 's/SHP/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
#sed -i 's/TBG/GLY/g' ${e3%.*}_r_prep1.pdb                 ###                    TBG GLY  T-BUTYL GLYCINE
#sed -i 's/3AH/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/DHI/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HIC/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HIS/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HIS
#sed -i 's/MHS/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/NEM/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/NEP/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HID/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  single delta N protonation
#sed -i 's/HIE/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                  single epsilon N protonation
#sed -i 's/3AH/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    3AH HIS
#sed -i 's/DDE/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    DDE HIS
#sed -i 's/DHI/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    DHI HIS  D-HISTIDINE
#sed -i 's/HIA/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HIA HIS  L-HISTIDINE AMIDE
#sed -i 's/HIC/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HIC HIS  4-METHYL-HISTIDINE
#sed -i 's/HIP/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HIP HIS  ND1-PHOSPHONOHISTIDINE...or commonly used doubly protonated state
#sed -i 's/HIQ/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HIQ HIS  MODIFIED HISTIDINE
#sed -i 's/HIS/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HIS HIS
#sed -i 's/HSO/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    HSO HIS  HISTIDINOL
#sed -i 's/MHS/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    MHS HIS  1-N-METHYLHISTIDINE
#sed -i 's/NEP/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    NEP HIS  N1-PHOSPHONOHISTIDINE
#sed -i 's/NZH/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    NZH HIS
#sed -i 's/OHI/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
#sed -i 's/PSH/HIS/g' ${e3%.*}_r_prep1.pdb                 ###                    PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
#sed -i 's/DIL/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
#sed -i 's/IIL/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
#sed -i 's/ILE/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                  ILE
#sed -i 's/B2I/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    B2I ILE  ISOLEUCINE BORONIC ACID
#sed -i 's/DIL/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    DIL ILE  D-ISOLEUCINE
#sed -i 's/IIL/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    IIL ILE  ISO-ISOLEUCINE
#sed -i 's/ILE/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    ILE ILE
#sed -i 's/ILX/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    ILX ILE  45-DIHYDROXYISOLEUCINE
#sed -i 's/IML/ILE/g' ${e3%.*}_r_prep1.pdb                 ###                    IML ILE  N-METHYLATED
#sed -i 's/ALY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/DLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/KCX/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LLP/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LYM/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LYS/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  LYS
#sed -i 's/LYZ/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/MLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/SHR/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/TRG/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/6CL/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    6CL LYS  6-CARBOXYLYSINE
#sed -i 's/ALY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    ALY LYS  N(6)-ACETYLLYSINE
#sed -i 's/API/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    API LYS  26-DIAMINOPIMELIC ACID
#sed -i 's/APK/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    APK LYS
#sed -i 's/AZK/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
#sed -i 's/B3K/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    B3K LYS  (3S)-37-DIAMINOHEPTANOIC ACID
#sed -i 's/BLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    BLY LYS  LYSINE BORONIC ACID
#sed -i 's/C1X/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    C1X LYS  MODIFIED LYSINE
#sed -i 's/CLG/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CLG LYS
#sed -i 's/CLH/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CLH LYS
#sed -i 's/CYJ/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    CYJ LYS  MODIFIED LYSINE
#sed -i 's/DLS/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    DLS LYS  DI-ACETYL-LYSINE
#sed -i 's/DLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    DLY LYS  D-LYSINE
#sed -i 's/DNL/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    DNL LYS  6-AMINO-HEXANAL
#sed -i 's/FHL/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    FHL LYS  MODIFIED LYSINE
#sed -i 's/GPL/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    GPL LYS  LYSINE GUANOSINE-5-MONOPHOSPHATE
#sed -i 's/IT1/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    IT1 LYS
#sed -i 's/KCX/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    KCX LYS  CARBAMOYLATED LYSINE
#sed -i 's/KGC/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    KGC LYS
#sed -i 's/KST/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
#sed -i 's/LA2/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LA2 LYS
#sed -i 's/LCK/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LCK LYS
#sed -i 's/LCX/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LCX LYS  CARBAMYLATED LYSINE
#sed -i 's/LDH/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LDH LYS  N~6~-ETHYL-L-LYSINE
#sed -i 's/LET/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LET LYS  ODIFIED LYSINE
#sed -i 's/LLP/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LLP LYS
#sed -i 's/LLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
#sed -i 's/LSO/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LSO LYS  MODIFIED LYSINE
#sed -i 's/LYM/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYM LYS  DEOXY-METHYL-LYSINE
#sed -i 's/LYN/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYN LYS  26-DIAMINO-HEXANOIC ACID AMIDE
#sed -i 's/LYP/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
#sed -i 's/LYR/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYR LYS  MODIFIED LYSINE
#sed -i 's/LYS/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYS LYS
#sed -i 's/LYX/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYX LYS  N-(2-COENZYME A)-PROPANOYL-LYSINE
#sed -i 's/LYZ/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    LYZ LYS  5-HYDROXYLYSINE
#sed -i 's/M2L/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    M2L LYS
#sed -i 's/M3L/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    M3L LYS  N-TRIMETHYLLYSINE
#sed -i 's/MCL/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
#sed -i 's/MLY/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    MLY LYS  METHYLATED LYSINE
#sed -i 's/MLZ/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    MLZ LYS  N-METHYL-LYSINE
#sed -i 's/OBS/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    OBS LYS  MODIFIED LYSINE
#sed -i 's/SLZ/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    SLZ LYS  L-THIALYSINE
#sed -i 's/XX1/LYS/g' ${e3%.*}_r_prep1.pdb                 ###                    XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE
#sed -i 's/BUG/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/CLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/DLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/LEU/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  LEU
#sed -i 's/MLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLN/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLP/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/1LU/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
#sed -i 's/2ML/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    2ML LEU  2-METHYLLEUCINE
#sed -i 's/BLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    BLE LEU  LEUCINE BORONIC ACID
#sed -i 's/BUG/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    BUG LEU  TERT-LEUCYL AMINE
#sed -i 's/CLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    CLE LEU  LEUCINE AMIDE
#sed -i 's/DCL/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
#sed -i 's/DLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    DLE LEU  D-LEUCINE
#sed -i 's/DNE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    DNE LEU  D-NORLEUCINE
#sed -i 's/DNG/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    DNG LEU  N-FORMYL-D-NORLEUCINE
#sed -i 's/DNM/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    DNM LEU  D-N-METHYL NORLEUCINE
#sed -i 's/FLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    FLE LEU  FUROYL-LEUCINE
#sed -i 's/HLU/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    HLU LEU  BETA-HYDROXYLEUCINE
#sed -i 's/LED/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    LED LEU  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/LEF/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    LEF LEU  2-5-FLUOROLEUCINE
#sed -i 's/LEU/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    LEU LEU
#sed -i 's/LNT/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    LNT LEU
#sed -i 's/MHL/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    MHL LEU  N-METHYLATED HYDROXY
#sed -i 's/MLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    MLE LEU  N-METHYLATED
#sed -i 's/MLL/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    MLL LEU  METHYL L-LEUCINATE
#sed -i 's/MNL/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    MNL LEU  4N-DIMETHYLNORLEUCINE
#sed -i 's/NLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    NLE LEU  NORLEUCINE
#sed -i 's/NLN/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    NLN LEU  NORLEUCINE AMIDE
#sed -i 's/NLO/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    NLO LEU  O-METHYL-L-NORLEUCINE
#sed -i 's/PLE/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    PLE LEU  LEUCINE PHOSPHINIC ACID
#sed -i 's/PPH/LEU/g' ${e3%.*}_r_prep1.pdb                 ###                    PPH LEU  PHENYLALANINE PHOSPHINIC ACID
#sed -i 's/CXM/MET/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/FME/MET/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/MET/MET/g' ${e3%.*}_r_prep1.pdb                 ###                  MET
#sed -i 's/MSE/MET/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/OMT/MET/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/AME/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    AME MET  ACETYLATED METHIONINE
#sed -i 's/CXM/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    CXM MET  N-CARBOXYMETHIONINE
#sed -i 's/ESC/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
#sed -i 's/FME/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    FME MET  FORMYL-METHIONINE
#sed -i 's/FOR/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    FOR MET
#sed -i 's/MET/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    MET MET
#sed -i 's/MHO/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    MHO MET  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/MME/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    MME MET  N-METHYL METHIONINE
#sed -i 's/MSE/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    MSE MET  ELENOMETHIONINE
#sed -i 's/MSO/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    MSO MET  METHIONINE SULFOXIDE
#sed -i 's/OMT/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    OMT MET  METHIONINE SULFONE
#sed -i 's/SME/MET/g' ${e3%.*}_r_prep1.pdb                 ###                    SME MET  METHIONINE SULFOXIDE
#sed -i 's/ASN/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                  ASN
#sed -i 's/MEN/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASN
#sed -i 's/AFA/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    AFA ASN  N-[7-METHYL-OCT-24-DIENOYL]ASPARAGINE
#sed -i 's/AHB/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    AHB ASN  BETA-HYDROXYASPARAGINE
#sed -i 's/ASN/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    ASN ASN
#sed -i 's/B3X/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    B3X ASN  (3S)-35-DIAMINO-5-OXOPENTANOIC ACID
#sed -i 's/DMH/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    DMH ASN  N4N4-DIMETHYL-ASPARAGINE
#sed -i 's/DSG/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    DSG ASN  D-ASPARAGINE
#sed -i 's/MEN/ASN/g' ${e3%.*}_r_prep1.pdb                 ###                    MEN ASN  GAMMA METHYL ASPARAGINE
#sed -i 's/DPR/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PRO
#sed -i 's/PRO/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                  PRO
#sed -i 's/1AB/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    1AB PRO  14-DIDEOXY-14-IMINO-D-ARABINITOL
#sed -i 's/2MT/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    2MT PRO
#sed -i 's/4FB/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    4FB PRO  (4S)-4-FLUORO-L-PROLINE
#sed -i 's/DPL/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    DPL PRO  4-OXOPROLINE
#sed -i 's/DPR/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    DPR PRO  D-PROLINE
#sed -i 's/H5M/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
#sed -i 's/HY3/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    HY3 PRO  3-HYDROXYPROLINE
#sed -i 's/HYP/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    HYP PRO  4-HYDROXYPROLINE
#sed -i 's/LPD/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    LPD PRO  L-PROLINAMIDE
#sed -i 's/P2Y/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
#sed -i 's/PCA/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    PCA PRO  5-OXOPROLINE
#sed -i 's/POM/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    POM PRO  CIS-5-METHYL-4-OXOPROLINE
#sed -i 's/PRO/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    PRO PRO
#sed -i 's/PRS/PRO/g' ${e3%.*}_r_prep1.pdb                 ###                    PRS PRO  THIOPROLINE
#sed -i 's/DGN/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLN
#sed -i 's/GLN/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                  GLN
#sed -i 's/DGN/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                    DGN GLN  D-GLUTAMINE
#sed -i 's/GHG/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                    GHG GLN  GAMMA-HYDROXY-GLUTAMINE
#sed -i 's/GLH/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                    GLH GLN
#sed -i 's/GLN/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                    GLN GLN
#sed -i 's/MGN/GLN/g' ${e3%.*}_r_prep1.pdb                 ###                    MGN GLN  2-METHYL-GLUTAMINE
#sed -i 's/ACL/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/AGM/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/ARG/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  ARG
#sed -i 's/ARM/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/DAR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/HAR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/HMR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/2MR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    2MR ARG  N3 N4-DIMETHYLARGININE
#sed -i 's/AAR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    AAR ARG  ARGININEAMIDE
#sed -i 's/ACL/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
#sed -i 's/AGM/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    AGM ARG  4-METHYL-ARGININE
#sed -i 's/ALG/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ALG ARG  GUANIDINOBUTYRYL GROUP
#sed -i 's/AR2/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
#sed -i 's/ARG/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ARG ARG
#sed -i 's/ARM/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ARM ARG  DEOXY-METHYL-ARGININE
#sed -i 's/ARO/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ARO ARG  C-GAMMA-HYDROXY ARGININE
#sed -i 's/BOR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    BOR ARG
#sed -i 's/CIR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    CIR ARG  CITRULLINE
#sed -i 's/DA2/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    DA2 ARG  MODIFIED ARGININE
#sed -i 's/DAR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    DAR ARG  D-ARGININE
#sed -i 's/HMR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    HMR ARG  BETA-HOMOARGININE
#sed -i 's/HRG/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    HRG ARG  L-HOMOARGININE
#sed -i 's/MAI/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    MAI ARG  DEOXO-METHYLARGININE
#sed -i 's/MGG/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    MGG ARG  MODIFIED D-ARGININE
#sed -i 's/NMM/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    NMM ARG  MODIFIED ARGININE
#sed -i 's/OPR/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    OPR ARG  C-(3-OXOPROPYL)ARGININE
#sed -i 's/ORQ/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    ORQ ARG  N~5~-ACETYL-L-ORNITHINE
#sed -i 's/TYZ/ARG/g' ${e3%.*}_r_prep1.pdb                 ###                    TYZ ARG  PARA ACETAMIDO BENZOIC ACID
#sed -i 's/DSN/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/MIS/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/OAS/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SAC/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SEL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SEP/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SER/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  SER
#sed -i 's/SET/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SVA/SER/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/B3S/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
#sed -i 's/BG1/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    BG1 SER
#sed -i 's/DHL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    DHL SER  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/DSE/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    DSE SER  D-SERINE N-METHYLATED
#sed -i 's/DSN/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    DSN SER  D-SERINE
#sed -i 's/FGP/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    FGP SER
#sed -i 's/GVL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
#sed -i 's/HSE/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    HSE SER  L-HOMOSERINE
#sed -i 's/HSL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    HSL SER  HOMOSERINE LACTONE
#sed -i 's/MC1/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    MC1 SER  METHICILLIN ACYL-SERINE
#sed -i 's/MIS/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    MIS SER  MODIFIED SERINE
#sed -i 's/N10/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
#sed -i 's/NC1/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    NC1 SER  NITROCEFIN ACYL-SERINE
#sed -i 's/OAS/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    OAS SER  O-ACETYLSERINE
#sed -i 's/OSE/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    OSE SER  O-SULFO-L-SERINE
#sed -i 's/PG1/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
#sed -i 's/PYR/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    PYR SER  CHEMICALLY MODIFIED
#sed -i 's/S1H/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
#sed -i 's/SAC/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SAC SER  N-ACETYL-SERINE
#sed -i 's/SBD/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SBD SER
#sed -i 's/SBG/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SBG SER  MODIFIED SERINE
#sed -i 's/SBL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SBL SER
#sed -i 's/SDP/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SDP SER
#sed -i 's/SEB/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SEB SER  O-BENZYLSULFONYL-SERINE
#sed -i 's/SEL/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SEL SER  2-AMINO-13-PROPANEDIOL
#sed -i 's/SEP/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SEP SER  E PHOSPHOSERINE
#sed -i 's/SER/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SER SER
#sed -i 's/SET/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SET SER  AMINOSERINE
#sed -i 's/SGB/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SGB SER  MODIFIED SERINE
#sed -i 's/SGR/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SGR SER  MODIFIED SERINE
#sed -i 's/SOY/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SOY SER  OXACILLOYL-ACYLATED SERINE
#sed -i 's/SUN/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SUN SER  TABUN CONJUGATED SERINE
#sed -i 's/SVA/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SVA SER  SERINE VANADATE
#sed -i 's/SVV/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SVV SER  MODIFIED SERINE
#sed -i 's/SVX/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SVX SER  MODIFIED SERINE
#sed -i 's/SVY/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SVY SER  MODIFIED SERINE
#sed -i 's/SVZ/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SVZ SER  MODIFIED SERINE
#sed -i 's/SXE/SER/g' ${e3%.*}_r_prep1.pdb                 ###                    SXE SER  MODIFIED SERINE
#sed -i 's/ALO/THR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/BMT/THR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/DTH/THR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/THR/THR/g' ${e3%.*}_r_prep1.pdb                 ###                  THR
#sed -i 's/TPO/THR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/AEI/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    AEI THR  ACYLATED THR
#sed -i 's/ALO/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    ALO THR  ALLO-THREONINE
#sed -i 's/BMT/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    BMT THR
#sed -i 's/CRO/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    CRO THR  CYCLIZED
#sed -i 's/CTH/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    CTH THR  4-CHLOROTHREONINE
#sed -i 's/DTH/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    DTH THR  D-THREONINE
#sed -i 's/OLT/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    OLT THR  O-METHYL-L-THREONINE
#sed -i 's/TBM/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    TBM THR
#sed -i 's/TH5/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    TH5 THR  O-ACETYL-L-THREONINE
#sed -i 's/THC/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    THC THR  N-METHYLCARBONYLTHREONINE
#sed -i 's/THR/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    THR THR
#sed -i 's/TMD/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    TMD THR  N-METHYLATED EPSILON C ALKYLATED
#sed -i 's/TPO/THR/g' ${e3%.*}_r_prep1.pdb                 ###                    TPO THR  HOSPHOTHREONINE
#sed -i 's/DIV/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/DVA/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/MVA/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/VAL/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                  VAL
#sed -i 's/B2V/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    B2V VAL  VALINE BORONIC ACID
#sed -i 's/DIV/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    DIV VAL  D-ISOVALINE
#sed -i 's/DVA/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    DVA VAL  D-VALINE
#sed -i 's/MNV/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    MNV VAL  N-METHYL-C-AMINO VALINE
#sed -i 's/MVA/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    MVA VAL  N-METHYLATED
#sed -i 's/NVA/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    NVA VAL  NORVALINE
#sed -i 's/VAD/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VAD VAL  DEAMINOHYDROXYVALINE
#sed -i 's/VAF/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VAF VAL  METHYLVALINE
#sed -i 's/VAL/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VAL VAL
#sed -i 's/VDL/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VDL VAL  (2R3R)-23-DIAMINOBUTANOIC ACID
#sed -i 's/VLL/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VLL VAL  (2S)-23-DIAMINOBUTANOIC ACID
#sed -i 's/VME/VAL/g' ${e3%.*}_r_prep1.pdb                 ###                    VME VAL  O- METHYLVALINE
#sed -i 's/DTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/HTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/LTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TPL/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TRO/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TRP/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                  TRP
#sed -i 's/BTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    BTR TRP  6-BROMO-TRYPTOPHAN
#sed -i 's/1TQ/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
#sed -i 's/23S/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    23S TRP  MODIFIED TRYPTOPHAN
#sed -i 's/32S/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    32S TRP  MODIFIED TRYPTOPHAN
#sed -i 's/32T/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    32T TRP  MODIFIED TRYPTOPHAN
#sed -i 's/4DP/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    4DP TRP
#sed -i 's/4FW/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    4FW TRP  4-FLUOROTRYPTOPHANE
#sed -i 's/4HT/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    4HT TRP  4-HYDROXYTRYPTOPHAN
#sed -i 's/4IN/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    4IN TRP  4-AMINO-L-TRYPTOPHAN
#sed -i 's/6CW/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    6CW TRP  6-CHLORO-L-TRYPTOPHAN
#sed -i 's/DTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    DTR TRP  D-TRYPTOPHAN
#sed -i 's/FTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    FTR TRP  FLUOROTRYPTOPHANE
#sed -i 's/HTR/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    HTR TRP  BETA-HYDROXYTRYPTOPHANE
#sed -i 's/PAT/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
#sed -i 's/TOX/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TOX TRP
#sed -i 's/TPL/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TPL TRP  TRYTOPHANOL
#sed -i 's/TQQ/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TQQ TRP
#sed -i 's/TRF/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRF TRP  N1-FORMYL-TRYPTOPHAN
#sed -i 's/TRN/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRN TRP  AZA-TRYPTOPHAN
#sed -i 's/TRO/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRO TRP  2-HYDROXY-TRYPTOPHAN
#sed -i 's/TRP/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRP TRP
#sed -i 's/TRQ/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRQ TRP
#sed -i 's/TRW/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRW TRP
#sed -i 's/TRX/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TRX TRP  6-HYDROXYTRYPTOPHAN
#sed -i 's/TTQ/TRP/g' ${e3%.*}_r_prep1.pdb                 ###                    TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
#sed -i 's/DTY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/IYR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/PAQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/PTR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/STY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYB/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYS/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  TYR
#sed -i 's/TYY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/1TY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    1TY TYR
#sed -i 's/2TY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    2TY TYR
#sed -i 's/3TY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    3TY TYR  MODIFIED TYROSINE
#sed -i 's/B3Y/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    B3Y TYR
#sed -i 's/CRQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    CRQ TYR
#sed -i 's/DBY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    DBY TYR  35 DIBROMOTYROSINE
#sed -i 's/DPQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    DPQ TYR  TYROSINE DERIVATIVE
#sed -i 's/DTY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    DTY TYR  D-TYROSINE
#sed -i 's/ESB/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    ESB TYR
#sed -i 's/FLT/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    FLT TYR  FLUOROMALONYL TYROSINE
#sed -i 's/FTY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
#sed -i 's/IYR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    IYR TYR  3-IODO-TYROSINE
#sed -i 's/MBQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    MBQ TYR
#sed -i 's/NIY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    NIY TYR  META-NITRO-TYROSINE
#sed -i 's/NBQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    NBQ TYR
#sed -i 's/OTY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    OTY TYR
#sed -i 's/PAQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    PAQ TYR  SEE REMARK 999
#sed -i 's/PTH/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
#sed -i 's/PTM/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
#sed -i 's/PTR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    PTR TYR  O-PHOSPHOTYROSINE
#sed -i 's/TCQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TCQ TYR  MODIFIED TYROSINE
#sed -i 's/TTS/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TTS TYR
#sed -i 's/TY2/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TY2 TYR  3-AMINO-L-TYROSINE
#sed -i 's/TY3/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TY3 TYR  3-HYDROXY-L-TYROSINE
#sed -i 's/TYB/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYB TYR  TYROSINAL
#sed -i 's/TYC/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYC TYR  L-TYROSINAMIDE
#sed -i 's/TYI/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYI TYR  35-DIIODOTYROSINE
#sed -i 's/TYN/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYN TYR  ADDUCT AT HYDROXY GROUP
#sed -i 's/TYO/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYO TYR
#sed -i 's/TYQ/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
#sed -i 's/TYR/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYR TYR
#sed -i 's/TYS/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYS TYR  INE SULPHONATED TYROSINE
#sed -i 's/TYT/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYT TYR
#sed -i 's/TYY/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
#sed -i 's/YOF/TYR/g' ${e3%.*}_r_prep1.pdb                 ###                    YOF TYR  3-FLUOROTYROSINE

#############################################################################################################################################################


pdb4amber -i ${poi%.*}_r_prep1.pdb -o ${poi%.*}_r_prep11.pdb
pdb4amber -i ${e3%.*}_r_prep1.pdb -o ${e3%.*}_r_prep11.pdb



sed -i '/UNL/!d' ${poi%.*}_ligand.pdb

sed -i '/UNL/!d' ${e3%.*}_ligand.pdb


sed -i '/HETATM.*ZN\|ZN.*HETATM/d' ${poi%.*}_ligand.pdb
sed -i '/HETATM.*ZN\|ZN.*HETATM/d' ${e3%.*}_ligand.pdb
sed -i '/HETATM.*NI\|NI.*HETATM/d' ${poi%.*}_ligand.pdb
sed -i '/HETATM.*NI\|NI.*HETATM/d' ${e3%.*}_ligand.pdb
sed -i '/HETATM.*MG\|MG.*HETATM/d' ${poi%.*}_ligand.pdb
sed -i '/HETATM.*MG\|MG.*HETATM/d' ${e3%.*}_ligand.pdb
sed -i '/HETATM.*CA\|CA.*HETATM/d' ${poi%.*}_ligand.pdb
sed -i '/HETATM.*CA\|CA.*HETATM/d' ${e3%.*}_ligand.pdb

sed -i 's/ATOM  /HETATM/g' ${poi%.*}_ligand.pdb
sed -i 's/ATOM  /HETATM/g' ${e3%.*}_ligand.pdb
sed -i 's/CL/Cl/g' ${poi%.*}_ligand.pdb
sed -i 's/CL/Cl/g' ${e3%.*}_ligand.pdb
sed -i 's/BR/Br/g' ${poi%.*}_ligand.pdb
sed -i 's/BR/Br/g' ${e3%.*}_ligand.pdb


#modified out oprot
#pdb4amber -i ${poi%.*}_ligand.pdb -o ${poi%.*}_ligand_oprot_prep1.pdb --reduce
#pdb4amber -i ${e3%.*}_ligand.pdb -o ${e3%.*}_ligand_oprot_prep1.pdb --reduce

#### EDITED
antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo mol2 -o ${poi%.*}_ligand_oprot_prep1_gas.mol2 -c gas -pf y

antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -dr n
#antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc 0 -pf y
else
	:
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -1 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +1 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -2 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +2 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -3 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +3 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +4 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -4 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -5 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +5 -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
	#pdb4amber -i ${poi%.*}_ligand.pdb -o ${poi%.*}_ligand_oprot_prep1.pdb --reduce
	#antechamber -fi pdb -i ${poi%.*}_ligand_oprot_prep1.pdb -fo mol2 -o ${poi%.*}_ligand_oprot_prep1_gas.mol2 -c gas -pf y
	antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -dr n
else
	:
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc 0 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -1 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +1 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -2 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +2 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -3 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +3 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +4 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -4 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -5 -pf y -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${poi%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +5 -pf y -dr n
else
        :
fi


if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas.mol2 ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas.mol2 ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 0
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 1
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -1
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 2
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -2
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 3
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -3
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas.mol2 ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas.mol2 ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 0 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 1 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -1 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 2 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -2 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 3 -dr n
else
        :
fi

if [ ! -f ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${poi%.*}_ligand.pdb -fo prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -3 -dr n
else
        :
fi

antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo mol2 -o ${e3%.*}_ligand_oprot_prep1_gas.mol2 -c gas -pf y
#antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -dr n
antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc 0 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -1 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +1 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -2 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +2 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -3 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +3 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +4 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -4 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -5 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +5 -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
	#grep "^CONECT" $e3 > ${e3%.*}_ligand2.pdb
	#cat ${e3%.*}_ligand2.pdb >> ${e3%.*}_ligand.pdb
        #pdb4amber -i ${e3%.*}_ligand.pdb -o ${e3%.*}_ligand_oprot_prep1.pdb --reduce
        #antechamber -fi pdb -i ${e3%.*}_ligand_oprot_prep1.pdb -fo mol2 -o ${e3%.*}_ligand_oprot_prep1_gas.mol2 -c gas -pf y
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc 0 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -1 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +1 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -2 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +2 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -3 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +3 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +4 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -4 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -5 -pf y -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${e3%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +5 -pf y -dr n
else
        :
fi


if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas.mol2 ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas.mol2 ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 0
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 1
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -1
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 2
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -2
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 3
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -3
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas.mol2 ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas.mol2 ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 0 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 1 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -1 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 2 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -2 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 3 -dr n
else
        :
fi

if [ ! -f ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${e3%.*}_ligand.pdb -fo prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -3 -dr n
else
        :
fi

#sed -i '/X.*nan\|nan.*X/d' ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi
#sed -i '/X.*nan\|nan.*X/d' ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi

parmchk2 -f prepi -i ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi -o ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
parmchk2 -f prepi -i ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi -o ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep11.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@



tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m2 = loadpdb ${e3%.*}_r_prep11.pdb        # load receptor with ion and cofactor
check m2
charge m2
saveamberparm m2 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

cp ${poi%.*}_r_prep11.pdb ${poi%.*}_r_prep00012.pdb
cp ${e3%.*}_r_prep11.pdb ${e3%.*}_r_prep00012.pdb


declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis000${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis000${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis000${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep009.pdb >> ${poi%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${poi%.*}_r_prep000${pcount}.pdb >> ${poi%.*}_r_prep000${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep000${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Parameter expansion substring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis000${count}.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis000${count}.out);
C="$(cut -d' ' -f1 <<<$del2)"
c=${C:1:3}
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep009.pdb >> ${e3%.*}_r_prep010.pdb
sed -e "/${c}.*${D}.*H/d" ${e3%.*}_r_prep000${pcount}.pdb >> ${e3%.*}_r_prep000${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep000${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

cp ${poi%.*}_r_prep00017.pdb ${poi%.*}_r_prep1112.pdb
cp ${e3%.*}_r_prep00017.pdb ${e3%.*}_r_prep1112.pdb

declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis11${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis11${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis11${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

sed -e "/${A}.*${B}.*H/d" ${poi%.*}_r_prep11${pcount}.pdb >> ${poi%.*}_r_prep11${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep11${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :

fi;
done;

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis11${count}.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis11${count}.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep006.pdb >> ${e3%.*}_r_prep007.pdb
sed -e "/${C}.*${D}.*H/d" ${e3%.*}_r_prep11${pcount}.pdb >> ${e3%.*}_r_prep11${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep11${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :

fi;
done;

cp ${poi%.*}_r_prep1145.pdb ${poi%.*}_r_prep11112.pdb
cp ${e3%.*}_r_prep1145.pdb ${e3%.*}_r_prep11112.pdb

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis111${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis111${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis111${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep009.pdb >> ${poi%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${poi%.*}_r_prep111${pcount}.pdb >> ${poi%.*}_r_prep111${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep111${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Parameter expansion substring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis111${count}.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis111${count}.out);
C="$(cut -d' ' -f1 <<<$del2)"
c=${C:1:3}
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep009.pdb >> ${e3%.*}_r_prep010.pdb
sed -e "/${c}.*${D}.*H/d" ${e3%.*}_r_prep111${pcount}.pdb >> ${e3%.*}_r_prep111${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep111${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

cp ${poi%.*}_r_prep11117.pdb ${poi%.*}_r_prep0012.pdb
cp ${e3%.*}_r_prep11117.pdb ${e3%.*}_r_prep0012.pdb

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do
if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis00${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis00${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis00${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep00${pcount}.pdb >> ${poi%.*}_r_prep00${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${poi%.*}_r_prep00${pcount}.pdb >> ${poi%.*}_r_prep00${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep00${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do
if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis00${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis00${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis00${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${e3%.*}_r_prep00${pcount}.pdb >> ${e3%.*}_r_prep00${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${e3%.*}_r_prep00${pcount}.pdb >> ${e3%.*}_r_prep00${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m2 = loadpdb ${e3%.*}_r_prep00${p2count}.pdb        # load receptor with ion and cofactor
check m2
charge m2
saveamberparm m2 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :

fi;
done;



if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O1P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis001.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis001.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis001.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep0047.pdb >> ${poi%.*}_r_prep001.pdb
sed -e "/O1P.*${A}.*${B}/d" ${poi%.*}_r_prep0047.pdb >> ${poi%.*}_r_prep001.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep001.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O1P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis002.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis002.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep0047.pdb >> ${e3%.*}_r_prep001.pdb
sed -e "/O1P.*${C}.*${D}/d" ${e3%.*}_r_prep0047.pdb >> ${e3%.*}_r_prep001.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep001.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O2P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis003.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis003.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis003.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep12.pdb >> ${poi%.*}_r_prep002.pdb
sed -e "/O2P.*${A}.*${B}/d" ${poi%.*}_r_prep001.pdb >> ${poi%.*}_r_prep002.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep002.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O2P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis004.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis004.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep12.pdb >> ${e3%.*}_r_prep002.pdb
sed -e "/O2P.*${C}.*${D}/d" ${e3%.*}_r_prep001.pdb >> ${e3%.*}_r_prep002.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep002.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O3P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis005.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis005.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis005.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep12.pdb >> ${poi%.*}_r_prep003.pdb
sed -e "/O3P.*${A}.*${B}/d" ${poi%.*}_r_prep002.pdb >> ${poi%.*}_r_prep003.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep003.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O3P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis006.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis006.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep12.pdb >> ${e3%.*}_r_prep003.pdb
sed -e "/O3P.*${C}.*${D}/d" ${e3%.*}_r_prep002.pdb >> ${e3%.*}_r_prep003.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep003.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;



if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis007.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis007.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis007.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep003.pdb >> ${poi%.*}_r_prep004.pdb
sed -e "/${A}.*${B}.*P/d" ${poi%.*}_r_prep003.pdb >> ${poi%.*}_r_prep004.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep004.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis008.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosi008.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep003.pdb >> ${e3%.*}_r_prep004.pdb
sed -e "/${C}.*${D}.*P/d" ${e3%.*}_r_prep003.pdb >> ${e3%.*}_r_prep004.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep004.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis009.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis009.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis009.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep004.pdb >> ${poi%.*}_r_prep005.pdb
sed -e "/${A}.*${B}.*P/d" ${poi%.*}_r_prep004.pdb >> ${poi%.*}_r_prep005.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep005.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis010.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis010.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep004.pdb >> ${e3%.*}_r_prep005.pdb
sed -e "/${C}.*${D}.*P/d" ${e3%.*}_r_prep004.pdb >> ${e3%.*}_r_prep005.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep005.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis011.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis011.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis011.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${poi%.*}_r_prep005.pdb >> ${poi%.*}_r_prep006.pdb
sed -e "/${A}.*${B}.*P/d" ${poi%.*}_r_prep005.pdb >> ${poi%.*}_r_prep006.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep006.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis012.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis012.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/P/d'  ${e3%.*}_r_prep005.pdb >> ${e3%.*}_r_prep006.pdb
sed -e "/${C}.*${D}.*P/d" ${e3%.*}_r_prep005.pdb >> ${e3%.*}_r_prep006.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep006.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis013.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis013.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis013.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep006.pdb >> ${poi%.*}_r_prep007.pdb
sed -e "/${A}.*${B}.*H/d" ${poi%.*}_r_prep006.pdb >> ${poi%.*}_r_prep007.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep007.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis014.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis014.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep006.pdb >> ${e3%.*}_r_prep007.pdb
sed -e "/${C}.*${D}.*H/d" ${e3%.*}_r_prep006.pdb >> ${e3%.*}_r_prep007.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep007.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis015.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis015.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis015.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep007.pdb >> ${poi%.*}_r_prep008.pdb
sed -e "/${A}.*${B}.*H/d" ${poi%.*}_r_prep007.pdb >> ${poi%.*}_r_prep008.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep008.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis016.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis016.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep007.pdb >> ${e3%.*}_r_prep008.pdb
sed -e "/${C}.*${D}.*H/d" ${e3%.*}_r_prep007.pdb >> ${e3%.*}_r_prep008.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep008.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis017.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis017.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis017.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep008.pdb >> ${poi%.*}_r_prep009.pdb
sed -e "/${A}.*${B}.*H/d" ${poi%.*}_r_prep008.pdb >> ${poi%.*}_r_prep009.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep009.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis018.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis018.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep008.pdb >> ${e3%.*}_r_prep009.pdb
sed -e "/${C}.*${D}.*H/d" ${e3%.*}_r_prep008.pdb >> ${e3%.*}_r_prep009.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep009.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis019.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis019.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis019.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep009.pdb >> ${poi%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${poi%.*}_r_prep009.pdb >> ${poi%.*}_r_prep010.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep010.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Parameter expansion substring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis020.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis020.out);
C="$(cut -d' ' -f1 <<<$del2)"
c=${C:1:3}
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep009.pdb >> ${e3%.*}_r_prep010.pdb
sed -e "/${c}.*${D}.*H/d" ${e3%.*}_r_prep009.pdb >> ${e3%.*}_r_prep010.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep010.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - HIS removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis021.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis021.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis021.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep010.pdb >> ${poi%.*}_r_prep011.pdb
sed -e "/HIS.*${B}.*H/d" ${poi%.*}_r_prep010.pdb >> ${poi%.*}_r_prep011.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep011.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - HIS removal in progress ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis022.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis022.out);
C="$(cut -d' ' -f1 <<<$del2)"
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep010.pdb >> ${e3%.*}_r_prep011.pdb
sed -e "/HIS.*${D}.*H/d" ${e3%.*}_r_prep010.pdb >> ${e3%.*}_r_prep011.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep011.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

declare -i count=23
declare -i pcount=11
declare -i p2count=12

for i in $(seq 35);
do
if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis0${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis0${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis0${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep0${pcount}.pdb >> ${poi%.*}_r_prep0${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${poi%.*}_r_prep0${pcount}.pdb >> ${poi%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
	:
fi;
done;

declare -i count=24
declare -i pcount=11
declare -i p2count=12

for i in $(seq 35);
do
if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis0${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis0${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis0${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${e3%.*}_r_prep0${pcount}.pdb >> ${e3%.*}_r_prep0${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${e3%.*}_r_prep0${pcount}.pdb >> ${e3%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m2 = loadpdb ${e3%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m2
charge m2
saveamberparm m2 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else 
	:

fi;
done;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P, O1P, O2P, O3P removal in progress ..."

sed -i '/ATOM.*\<P\>/d' ${e3%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O1P\>/d' ${e3%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O2P\>/d' ${e3%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O3P\>/d' ${e3%.*}_r_prep046.pdb



tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep046.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@
fi;

if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P, O1P, O2P, O3P removal in progress ..."

sed -i '/ATOM.*\<P\>/d' ${poi%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O1P\>/d' ${poi%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O2P\>/d' ${poi%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O3P\>/d' ${poi%.*}_r_prep046.pdb



tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep046.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@
fi;

declare -i count=1
declare -i pcount=46
declare -i p2count=47

for i in $(seq 8);
do
if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis22${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis22${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis22${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${poi%.*}_r_prep009.pdb >> ${poi%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${poi%.*}_r_prep0${pcount}.pdb >> ${poi%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${poi%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

declare -i count=1
declare -i pcount=46
declare -i p2count=47

for i in $(seq 8);
do
if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Parameter expansion substring extraction ..."
grep "^FATAL" ./slurm-${SLURM_JOB_ID}.out | tail -1 > diagnosis22${count}.out;
del2=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis22${count}.out);
C="$(cut -d' ' -f1 <<<$del2)"
c=${C:1:3}
D="$(cut -d' ' -f2 <<<$del2)"

#sed -e '/${C}/!b' -e '/${D}/!b' -e '/H/d'  ${e3%.*}_r_prep009.pdb >> ${e3%.*}_r_prep010.pdb
sed -e "/${c}.*${D}.*H/d" ${e3%.*}_r_prep0${pcount}.pdb >> ${e3%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;


if [ ! -f ${poi%.*}_r_prep2.inpcrd ] || [ ! -s ${poi%.*}_r_prep2.inpcrd ]; then

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${poi%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${poi%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${poi%.*}_r_special2.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${poi%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

fi;

if [ ! -f ${e3%.*}_r_prep2.inpcrd ] || [ ! -s ${e3%.*}_r_prep2.inpcrd ]; then

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${e3%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameter
loadamberparams ${e3%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${e3%.*}_r_special2.pdb        # load receptor with ion and cofactor
check m1
charge m1
saveamberparm m1 ${e3%.*}_r_prep2.prmtop ${e3%.*}_r_prep2.inpcrd
quit
@

fi;

ambpdb -p ${poi%.*}_r_prep2.prmtop -c ${poi%.*}_r_prep2.inpcrd -mol2 -sybyl > ${poi%.*}_r_prep3.mol2
ambpdb -p ${e3%.*}_r_prep2.prmtop -c ${e3%.*}_r_prep2.inpcrd -mol2 -sybyl > ${e3%.*}_r_prep3.mol2

echo "STARTING PROPOSE PP DOCKING"

#bash push_renum.sh ${poi%.*}_r_prep1.pdb ${poi%.*}_r_prep1_renum.txt ${poi%.*}_r_prep3.mol2 > ${poi%.*}_r_prep3_renum.mol2
#if [ ! -f ${poi%.*}_r_prep3_renum.mol2 ] || [ ! -s ${poi%.*}_r_prep3_renum.mol2 ]; then
#        echo "contournement de renum!!!"
#        awk '!($2="")' ${poi%.*}_r_prep1_renum.txt >> ${poi%.*}_f_prep1_renum.txt
#        bash push_renum.sh ${poi%.*}_r_prep1.pdb ${poi%.*}_f_prep1_renum.txt ${poi%.*}_r_prep3.mol2 > ${poi%.*}_r_prep3_renum.mol2
#fi
#
#if [ ! -f ${poi%.*}_r_prep3_renum.mol2 ] || [ ! -s ${poi%.*}_r_prep3_renum.mol2 ]; then
#        echo "error is persisting in push_renum.sh for $poi"
#fi

#bash push_renum.sh ${e3%.*}_r_prep1.pdb ${e3%.*}_r_prep1_renum.txt ${e3%.*}_r_prep3.mol2 > ${e3%.*}_r_prep3_renum.mol2
#if [ ! -f ${e3%.*}_r_prep3_renum.mol2 ] || [ ! -s ${e3%.*}_r_prep3_renum.mol2 ]; then
#        awk '!($2="")' ${e3%.*}_r_prep1_renum.txt >> ${e3%.*}_f_prep1_renum.txt
#        bash push_renum.sh ${e3%.*}_r_prep1.pdb ${poi%.*}_f_prep1_renum.txt ${e3%.*}_r_prep3.mol2 > ${e3%.*}_r_prep3_renum.mol2
#fi
#
#if [ ! -f ${e3%.*}_r_prep3_renum.mol2 ] || [ ! -s ${e3%.*}_r_prep3_renum.mol2 ]; then
#        echo "error is persisting in push_renum.sh for $e3"
#fi


sed -i '/HETATM.*ZN\|ZN.*HETATM/d' ${poi%.*}_r_prep11.pdb
sed -i '/HETATM.*ZN\|ZN.*HETATM/d' ${e3%.*}_r_prep11.pdb
sed -i '/HETATM.*NI\|NI.*HETATM/d' ${poi%.*}_r_prep11.pdb
sed -i '/HETATM.*NI\|NI.*HETATM/d' ${e3%.*}_r_prep11.pdb
sed -i '/HETATM.*MG\|MG.*HETATM/d' ${poi%.*}_r_prep11.pdb
sed -i '/HETATM.*MG\|MG.*HETATM/d' ${e3%.*}_r_prep11.pdb
sed -i '/HETATM.*CA\|CA.*HETATM/d' ${poi%.*}_r_prep11.pdb
sed -i '/HETATM.*CA\|CA.*HETATM/d' ${e3%.*}_r_prep11.pdb

python get_pocket.py ${poi%.*}_r_prep11.pdb ${poi%.*}_r_pocket.txt
python get_pocket.py ${e3%.*}_r_prep11.pdb ${e3%.*}_r_pocket.txt

sed -i '$ s/.$//' ${poi%.*}_r_pocket.txt
sed -i '$ s/.$//' ${e3%.*}_r_pocket.txt

sed -i '1s/^/:/' ${poi%.*}_r_pocket.txt
sed -i '1s/^/:/' ${e3%.*}_r_pocket.txt

pocket=$(head -n 1 ${poi%.*}_r_pocket.txt)
pocket2=$(head -n 1 ${e3%.*}_r_pocket.txt)

echo "$pocket"
echo "$pocket2"

ambmask -p ${poi%.*}_r_prep2.prmtop -c ${poi%.*}_r_prep2.inpcrd -prnlev 0 -out pdb -find "$pocket"| awk '/ATOM/{print $2}' > ${poi%.*}_r_prep3.hit
ambmask -p ${e3%.*}_r_prep2.prmtop -c ${e3%.*}_r_prep2.inpcrd -prnlev 0 -out pdb -find "$pocket2"| awk '/ATOM/{print $2}' > ${e3%.*}_r_prep3.hit


#BIN=${PROPOSEHOME?"undefined"}/bin

./ProPOSE TARGET=${poi%.*}_r_prep3.mol2 \
          LIGAND=${e3%.*}_r_prep3.mol2 \
          OUTTAR=${poi%.*}_r_prep_pp.mol2 \
          OUTLIG=${e3%.*}_l_prep_pp.mol2 \
          TARGET_PRM=${poi%.*}_r_prep2.prmtop \
          LIGAND_PRM=${e3%.*}_r_prep2.prmtop \
	  LIGAND_RMS=${e3%.*}_r_prep3.mol2 \
	  LICENSE=ProPOSE.lic \
          NOUT=100 \
	  NUM_THREAD=32
#	  
	
if [ -f ${poi%.*}_r_prep_pp.mol2 ] || [ -s ${poi%.*}_r_prep_pp.mol2 ]; then
	rm *special* *hit *prep2* *prep1* *ligand* sqm.* ANTECHAMBER* ATOMTYPE* *log *noHET* *sed* *sslink* *txt *renum*mol2
	echo "SUCCESS"
else
	echo "FAILED"
fi



rm *scores pp_scores.csv
awk '/SCORES:/,/Saved/' slurm*out  >> trial.scores
tail -n +2 trial.scores >> trial2.scores
head -n -1 trial2.scores > temp.txt ; mv temp.txt trial2.scores
column -t trial2.scores >> trial3.scores
sed -e 's/\s\+/,/g' trial3.scores > pp_scores.csv
#sed -i '2,102s/^.//' pp_scores.csv
rm *scores


rm *filtered* *fix*
for p in *pp*mol2; do sed -e '/P.*UNL/!b' -e 's/C.3/P.3/g' $p >> ${p%.*}_p1fixed.mol2;done

python filter_mol2.py


sbatch moe.sh $poi $e3

