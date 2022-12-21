#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=minimize1

source ~/.bashrc
conda activate amber

poi=$1
linker=$2
link=${linker##*/}

cat $poi $linker >> ${poi%.*}_${link%.*}_linker_cat.pdb
pdb4amber -i ${poi%.*}_${link%.*}_linker_cat.pdb -o amb_${poi%.*}_${link%.*}_linker_cat.pdb

python unclash.py amb_${poi%.*}_${link%.*}_linker_cat.pdb

cp unclashed_amb_${poi%.*}_${link%.*}_linker_cat.pdb amb_${poi%.*}_${link%.*}.pdb

mkdir amb_${poi%.*}_${link%.*}
mv amb_${poi%.*}_${link%.*}.pdb amb_${poi%.*}_${link%.*}
cd amb_${poi%.*}_${link%.*}
$SCHRODINGER/utilities/prepwizard -NOJOBID -LOCAL -noepik amb_${poi%.*}_${link%.*}.pdb mae_amb_${poi%.*}_${link%.*}_mae.pdb -f 3 -WAIT -j STAGE1_${poi%.*}_${link%.*} -rmsd 5.0

##ADD CONECT, REMOVE HYD99
#cp mae_amb_${poi%.*}_${link%.*}_mae.pdb inspect_${poi%.*}_${link%.*}.pdb

sed -i '/END/d' mae_amb_${poi%.*}_${link%.*}_mae.pdb;
sed -i '/MASTER/d' mae_amb_${poi%.*}_${link%.*}_mae.pdb;

grep -n -m 2 -e C90 -e O90 -e N90 -e P90 -e S90 mae_amb_${poi%.*}_${link%.*}_mae.pdb >> ${poi%.*}_${link%.*}_in.file
grep -n -m 2 -e C91 -e O91 -e N91 -e P91 -e S91 mae_amb_${poi%.*}_${link%.*}_mae.pdb >> ${poi%.*}_${link%.*}_in2.file
grep -n -m 2 -e C92 -e O92 -e N92 -e P92 -e S92 mae_amb_${poi%.*}_${link%.*}_mae.pdb >> ${poi%.*}_${link%.*}_in3.file
grep -n -m 2 -e C93 -e O93 -e N93 -e P93 -e S93 mae_amb_${poi%.*}_${link%.*}_mae.pdb >> ${poi%.*}_${link%.*}_in4.file
grep -n -m 2 -e C94 -e O94 -e N94 -e P94 -e S94 mae_amb_${poi%.*}_${link%.*}_mae.pdb >> ${poi%.*}_${link%.*}_in5.file
tr -d ' ' < ${poi%.*}_${link%.*}_in.file | cut -d':' -f1 | sort -u > ${poi%.*}_${link%.*}_out.file
tr -d ' ' < ${poi%.*}_${link%.*}_in2.file | cut -d':' -f1 | sort -u > ${poi%.*}_${link%.*}_out2.file
tr -d ' ' < ${poi%.*}_${link%.*}_in3.file | cut -d':' -f1 | sort -u > ${poi%.*}_${link%.*}_out3.file
tr -d ' ' < ${poi%.*}_${link%.*}_in4.file | cut -d':' -f1 | sort -u > ${poi%.*}_${link%.*}_out4.file
tr -d ' ' < ${poi%.*}_${link%.*}_in5.file | cut -d':' -f1 | sort -u > ${poi%.*}_${link%.*}_out5.file
var1=$(awk 'NR==1' ${poi%.*}_${link%.*}_out.file);
var2=$(awk 'NR==1' ${poi%.*}_${link%.*}_out2.file);
var3=$(awk 'NR==1' ${poi%.*}_${link%.*}_out3.file);
var4=$(awk 'NR==1' ${poi%.*}_${link%.*}_out4.file);
var5=$(awk 'NR==1' ${poi%.*}_${link%.*}_out5.file);
#var6=$(awk 'NR==1' ${poi%.*}_${link%.*}_out5.file);

var1=$((var1-3));
var2=$((var2-3));
var3=$((var3-3));
var4=$((var4-3));
var5=$((var5-3));

#var6=$((var6-3));

if [[ "$var1" -eq -3 ]] && [[ "$var5" -ne -3 ]];
then
	:
elif [[ "$var2" -ne -3 ]] && [[ "$var3" -eq -3 ]];
then
	if [ ${#var1} -gt 4 ] && [ ${#var2} -gt 4 ];
	then
		string="CONECT${var1}${var2}";
	elif [ ${#var1} -gt 4 ] && [ ${#var2} -lt 5 ];
	then 
		string="CONECT${var1} ${var2}";
	elif [ ${#var2} -gt 4 ] && [ ${#var1} -lt 5 ];
	then 
		string="CONECT ${var1}${var2}";
	else 
		string="CONECT ${var1} ${var2}";
	fi;
	sed -i '0,/H90.*UNL/{/H90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H901.*UNL/{/H901.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H91.*UNL/{/H91.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H911.*UNL/{/H911.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H92.*UNL/{/H92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/H921.*UNL/{/H921.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HO90.*UNL/{/HO90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HO91.*UNL/{/HO91.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HO92.*UNL/{/HO92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HN90.*UNL/{/HN90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HN91.*UNL/{/HN91.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HN92.*UNL/{/HN92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HS90.*UNL/{/HS90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HS91.*UNL/{/HS91.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HS92.*UNL/{/HS92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HP90.*UNL/{/HP90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HP91.*UNL/{/HP91.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HP92.*UNL/{/HP92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
elif [[ "$var3" -ne -3 ]] && [[ "$var4" -eq -3 ]];
then 
	if [ ${#var1} -gt 4 ] && [ ${#var2} -gt 4 ] && [ ${#var3} -gt 4 ];
        then
                string="CONECT${var1}${var2}${var3}";
        elif [ ${#var1} -lt 5 ] && [ ${#var2} -gt 4 ] && [ ${#var3} -gt 4 ];
        then
                string="CONECT ${var1}${var2}${var3}";
        elif [ ${#var1} -lt 5 ] && [ ${#var2} -lt 5 ] && [ ${#var3} -gt 4 ];
        then
                string="CONECT ${var1} ${var2}${var3}";
        elif [ ${#var1} -gt 4 ] && [ ${#var2} -lt 5 ] && [ ${#var3} -gt 4 ];
	then
                string="CONECT${var1} ${var2}${var3}";
	elif [ ${#var1} -gt 4 ] && [ ${#var2} -lt 5 ] && [ ${#var3} -lt 5 ];
	then
                string="CONECT${var1} ${var2} ${var3}";
	elif [ ${#var1} -lt 5 ] && [ ${#var2} -gt 4 ] && [ ${#var3} -lt 5 ];
        then
                string="CONECT ${var1}${var2} ${var3}";
	elif [ ${#var1} -gt 4 ] && [ ${#var2} -gt 4 ] && [ ${#var3} -lt 5 ];
        then
                string="CONECT${var1}${var2} ${var3}";
	else
		string="CONECT ${var1} ${var2} ${var3}";
	fi;
        sed -i '0,/H90.*UNL/{/H90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H901.*UNL/{/H901.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H902.*UNL/{/H902.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H91.*UNL/{/H91.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/H911.*UNL/{/H911.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H92.*UNL/{/H92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/H921.*UNL/{/H921.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HO90.*UNL/{/HO90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HO91.*UNL/{/HO91.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HO92.*UNL/{/HO92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HN90.*UNL/{/HN90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HN91.*UNL/{/HN91.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HN92.*UNL/{/HN92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HS^0.*UNL/{/HS90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HS91.*UNL/{/HS91.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HS92.*UNL/{/HS92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HP90.*UNL/{/HP90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HP91.*UNL/{/HP91.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HP92.*UNL/{/HP92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
elif [[ "$var4" -ne -3 ]];
then
        if [ ${#var1} -gt 4 ] && [ ${#var4} -gt 4 ] && [ ${#var3} -gt 4 ];
        then
                string="CONECT${var1}${var4}${var4}${var3}";
        elif [ ${#var1} -lt 5 ] && [ ${#var4} -gt 4 ] && [ ${#var3} -gt 4 ];
        then
                string="CONECT ${var1}${var4}${var4}${var3}";
        elif [ ${#var1} -lt 5 ] && [ ${#var4} -lt 5 ] && [ ${#var3} -gt 4 ];
        then
                string="CONECT ${var1} ${var4} ${var4}${var3}";
        elif [ ${#var1} -gt 4 ] && [ ${#var4} -lt 5 ] && [ ${#var3} -gt 4 ];
        then
                string="CONECT${var1} ${var4} ${var4}${var3}";
        elif [ ${#var1} -gt 4 ] && [ ${#var4} -lt 5 ] && [ ${#var3} -lt 5 ];
        then
                string="CONECT${var1} ${var4} ${var4} ${var3}";
        elif [ ${#var1} -lt 5 ] && [ ${#var4} -gt 4 ] && [ ${#var3} -lt 5 ];
        then
                string="CONECT ${var1}${var4}${var4} ${var3}";
        elif [ ${#var1} -gt 4 ] && [ ${#var4} -gt 4 ] && [ ${#var3} -lt 5 ];
        then
                string="CONECT${var1}${var4}${var4} ${var3}";
        else
                string="CONECT ${var1} ${var4} ${var4} ${var3}";
        fi;
	sed -i '0,/H90.*UNL/{/H90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H901.*UNL/{/H901.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H91.*UNL/{/H91.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H902.*UNL/{/H902.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H903.*UNL/{/H903.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/HO90.*UNL/{/HO90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/HN90.*UNL/{/HN90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/HN90.*UNL/{/HN90.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H93.*UNL/{/H93.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/H931.*UNL/{/H931.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/H932.*UNL/{/H932.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HO93.*UNL/{/HO93.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/HN93.*UNL/{/HN93.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HN93.*UNL/{/HN93.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
	sed -i '0,/H92.*UNL/{/H92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/H921.*UNL/{/H921.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HO92.*UNL/{/HO92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HN92.*UNL/{/HN92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HS92.*UNL/{/HS92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
        sed -i '0,/HP92.*UNL/{/HP92.*UNL/d;}' mae_amb_${poi%.*}_${link%.*}_mae.pdb
else
	echo "combination not found"
fi;
	
#string2="CONECT ${var3} ${var4}"

echo $string >> mae_amb_${poi%.*}_${link%.*}_mae.pdb;
#if [ -z "$string2" ]
#then
#        :
#else
#        echo $string2 >> mae_amb_${poi%.*}_${link%.*}_mae.pdb;
#fi
echo "END" >> mae_amb_${poi%.*}_${link%.*}_mae.pdb;
rm ${poi%.*}_${link%.*}*file

mkdir mae_amb_${poi%.*}_${link%.*}_mae
mv mae_amb_${poi%.*}_${link%.*}_mae.pdb mae_amb_${poi%.*}_${link%.*}_mae
rm *
cd mae_amb_${poi%.*}_${link%.*}_mae

$SCHRODINGER/utilities/prepwizard -NOJOBID -LOCAL -noepik mae_amb_${poi%.*}_${link%.*}_mae.pdb nhyd99_mae_amb_${poi%.*}_${link%.*}.pdb -noccd -f 3 -WAIT -j STAGE1_${poi%.*}_${link%.*} -rmsd 5.0

cp ../../utilities/valency0.py ./

if [ ! -f nhyd99_mae_amb_${poi%.*}_${link%.*}.pdb ]; then
python valency0.py mae_amb_${poi%.*}_${link%.*}_mae.pdb
$SCHRODINGER/utilities/prepwizard -NOJOBID -LOCAL -noepik valency_mae_amb_${poi%.*}_${link%.*}_mae.pdb nhyd99_mae_amb_${poi%.*}_${link%.*}.pdb -f 3 -WAIT -j STAGE1_${poi%.*}_${link%.*} -rmsd 5.0
fi

mkdir nhyd99_mae_amb_${poi%.*}_${link%.*}
mv nhyd99_mae_amb_${poi%.*}_${link%.*}.pdb nhyd99_mae_amb_${poi%.*}_${link%.*}
cd nhyd99_mae_amb_${poi%.*}_${link%.*}
$SCHRODINGER/utilities/prepwizard -NOJOBID -LOCAL -noepik nhyd99_mae_amb_${poi%.*}_${link%.*}.pdb out_nhyd99_mae_amb_${poi%.*}_${link%.*}.pdb -f 3 -WAIT -j STAGE1_${poi%.*}_${link%.*} -rmsd 5.0


sed -i '/HETATM/s/  BR  / BR   /g' out_nhyd99_mae_amb_${poi%.*}_${link%.*}.pdb
sed -i '/HETATM/s/  CL  / CL   /g' out_nhyd99_mae_amb_${poi%.*}_${link%.*}.pdb


#rm nhyd99_mae_amb_${poi%.*}_${link%.*}.pdb mae_amb_${poi%.*}_${link%.*}.pdb amb_${poi%.*}_${link%.*}.pdb *sslink* *nonprot* *renum*

mv out_nhyd99_mae_amb_${poi%.*}_${link%.*}.pdb ../../../
rm *
cd ../../../
grep "^HETATM" out_nhyd99_mae_amb_${poi%.*}_${link%.*}.pdb > ${poi%.*}_${link%.*}_ligand_comparison.pdb
grep "^HETATM" $poi > ${poi%.*}_${link%.*}_ref_ligand.pdb

obrms ${poi%.*}_${link%.*}_ref_ligand.pdb ${poi%.*}_${link%.*}_ligand_comparison.pdb >> ${poi%.*}_${link%.*}_rmsd.rms

#rm -r linker_library/temp


