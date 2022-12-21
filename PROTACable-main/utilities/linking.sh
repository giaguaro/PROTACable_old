#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=PROTACs

source ~/.bashrc
conda activate py39

rm *pp*pdb *sslink* *renum* *nonprot*
for m1 in poi*_prep_pp.mol2; do 
	obabel -imol2 $m1 -opdb -O ${m1%.*}.pdb -m
	for p1 in poi*_r_prep_pp*pdb; do
		sed -i -e '/UNL/ s/ATOM  /HETATM/' $p1;
		sed -i '/CONECT/d' $p1;
		sed -i '/MASTER/d' $p1;
		sed -i '/END/d' $p1;
		grep "^HETATM" $p1 > ${p1%.*}_ligand.pdb;
		sed -i '/HETATM/d' $p1;
	done;
done;
for m2 in e3*_prep_pp.mol2; do
        obabel -imol2 $m2 -opdb -O ${m2%.*}.pdb -m
	for p2 in e3*_l_prep_pp*pdb; do
                sed -i -e '/UNL/ s/ATOM  /HETATM/' $p2;
		sed -i '/CONECT/d' $p2;
                sed -i '/MASTER/d' $p2;
                sed -i '/END/d' $p2;
		grep "^HETATM" $p2 > ${p2%.*}_ligand.pdb;
		sed -i '/HETATM/d' $p2;
		v1=$(echo "$p2" | cut -d '_' -f 9);
		v1s=$(echo "$v1" | cut -d '.' -f 1);
		echo "$v1s";
		cat $p2 >> poi*pp*${v1}.pdb;
		cat poi*${v1s}_ligand.pdb >> poi*pp*${v1}.pdb;
		cat ${p2%.*}_ligand.pdb >> poi*pp*${v1}.pdb;
		rm $p2;
		mv poi*pp*${v1}.pdb poi_${p2%.*}.pdb;
		obabel -ipdb poi_${p2%.*}.pdb -opdb -O poi_${p2%.*}_complex.pdb;
                pdb4amber -i poi_${p2%.*}_complex.pdb -o poi_${p2%.*}_complex_clean.pdb;
                rm poi_${p2%.*}_complex.pdb poi_${p2%.*}.pdb *sslink* *renum* *nonprot* poi*${v1s}_ligand.pdb ${p2%.*}_ligand.pdb
		grep -n -m 2 -e C99 -e O99 -e N99 poi_${p2%.*}_complex_clean.pdb >> in.file
		tr -d ' ' < in.file | cut -d':' -f1 | sort -u > out.file
		var1=$(awk 'NR==1' out.file);
		var2=$(awk 'NR==2' out.file);
		var1=$((var1-1));
		var2=$((var2-1));
		sed -i '/END/d' poi_${p2%.*}_complex_clean.pdb;
		string="CONECT ${var1} ${var2}"
		echo $string >> poi_${p2%.*}_complex_clean.pdb;
		echo "END" >> poi_${p2%.*}_complex_clean.pdb;
		rm in.file out.file
	done;
done;
