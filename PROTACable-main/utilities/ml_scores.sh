#!/bin/bash

#SBATCH --partition=normal
#SBATCH --job-name=scores
#SBATCH --ntasks=10
source ~/.bashrc
conda activate py39

##input is for p in *xx01-out-out-1.pdb; do sbatch 
#new
#pv_out_nhyd99_mae_complex_1err_33_6GFZ_6_pp_model_xx14-out-out-1.pdb

#for p in *_xx01-out-out-1.pdb;
#do
in=$1
v2=$(echo "$in" | cut -d '_' -f 6-9)


for input in pv_out_nhyd99_mae_complex_${v2}_*-1.pdb;
do
v2=$(echo "$input" | cut -d '_' -f 6-9)
v3=$(echo "$input" | cut -d '_' -f 12 | cut -d '-' -f 1)

original=$PWD

#ligand_amide_mpro_docked.sdf
dir_lig_poi_dock="../"

#Cereblon_3QQU_1_docked2_0.sdf
dir_lig_e3_dock="../utilities/e3_scores_sdf"

#poi_eraf2_cp_rec_6CB5_1_linker_docked_short_rmsd.rms
dir_rmsd_linker="../"

#e3_vhl_5t8e_49_docked_rec_6GFZ_6_pair
dir_propose_dock="../"

#pv_out_nhyd99_mae_complex_1err_33_6GFZ_6_pp_model_xx15-out-out.csv
dir_mmgbsa_score="./"

##get POI dock scores

cp get_scores.py ${dir_lig_poi_dock}
cd ${dir_lig_poi_dock} 
#v1=$(echo "$input" | cut -d '_' -f 6,7 | tr '[:lower:]' '[:upper:]')
v1=$(ls *_docked.sdf)
python get_scores.py $v1 ${v1%.*}_scores.txt ${v1%.*}_scores.txt
icm_poi_lig_score=$(head -2 ${v1%.*}_scores.txt | tail -1 | cut -d ',' -f 2)

cd ${original}

##get E3 dock scores

cp get_scores.py ${dir_lig_e3_dock}
cd ${dir_lig_e3_dock}
v4=$(echo "$input" | cut -d '_' -f 6-7 | tr '[:lower:]' '[:upper:]')
python get_scores.py *${v4}_docked2_0.sdf e3_${v4}_scores.txt e3_${v4}_scores.txt
icm_e3_lig_score=$(head -2 e3_${v4}_scores.txt | tail -1 | cut -d ',' -f 2)

cd ${original}


##get rmsd
v5=$(echo "$input" | cut -d '_' -f 8-9)
cd ${dir_rmsd_linker}
link_rms=$(ls *_${v5}_linker_docked_short_rmsd.rms)
rmsd_linker=$(cat ${link_rms} | awk '{ print $3 }')

cd ${original}

##get propose scores

v6=${v3:2:4}
#v5=$(echo "$in" | cut -d '_' -f 12)
v7=$(echo "$input" | cut -d '_' -f 6-7)
cd ${dir_propose_dock}/e3_*_${v7}_docked_rec_${v5}_pair
#rank,score,csiz,vdw,ele,ehb,flw,bpa,bnpa,dint,dcost,erf,rms,Pose
propose_score=$(cat pp_scores_filtered.csv| awk -F, '{ print $2 }'|sed "$((10#${v6}+1))q;d")
propose_csiz=$(cat pp_scores_filtered.csv| awk -F, '{ print $3 }'|sed "$((10#${v6}+1))q;d")
propose_vdw=$(cat pp_scores_filtered.csv| awk -F, '{ print $4 }'|sed "$((10#${v6}+1))q;d")
propose_ele=$(cat pp_scores_filtered.csv| awk -F, '{ print $5 }'|sed "$((10#${v6}+1))q;d")
propose_ehb=$(cat pp_scores_filtered.csv| awk -F, '{ print $6 }'|sed "$((10#${v6}+1))q;d")
propose_flw=$(cat pp_scores_filtered.csv| awk -F, '{ print $7 }'|sed "$((10#${v6}+1))q;d")
propose_bpa=$(cat pp_scores_filtered.csv| awk -F, '{ print $8 }'|sed "$((10#${v6}+1))q;d")
propose_bnpa=$(cat pp_scores_filtered.csv| awk -F, '{ print $9 }'|sed "$((10#${v6}+1))q;d")
propose_dint=$(cat pp_scores_filtered.csv| awk -F, '{ print $10 }'|sed "$((10#${v6}+1))q;d")
propose_dcost=$(cat pp_scores_filtered.csv| awk -F, '{ print $11 }'|sed "$((10#${v6}+1))q;d")
propose_erf=$(cat pp_scores_filtered.csv| awk -F, '{ print $12 }'|sed "$((10#${v6}+1))q;d")

cd ${original}


##get mmgbsa scores
cp scores_ml.py ${dir_mmgbsa_score}

cd ${dir_mmgbsa_score}

echo "${v4}_${v5}_${v6} $icm_e3_lig_score $icm_e3_lig_score $rmsd_linker $propose_score $propose_csiz $propose_vdw $propose_ele $propose_ehb $propose_flw $propose_bpa $propose_bnpa $propose_dint $propose_dcost $propose_erf pv_out_nhyd99_mae_complex_${v2}_pp_model_xx${v6}-out-out.csv"

python scores_ml.py ${v2}_${v6} $icm_poi_lig_score $icm_e3_lig_score $rmsd_linker $propose_score $propose_csiz $propose_vdw $propose_ele $propose_ehb $propose_flw $propose_bpa $propose_bnpa $propose_dint $propose_dcost $propose_erf pv_out_nhyd99_mae_complex_${v2}_pp_model_xx${v6}-out-out.csv


done;


python combine_ml_scores.py ${v2}

mv ${v2}_all_ml_scores.csv ../ml_scores

rm ${v2}_*ml_scores.csv

cd ../ml_scores

python final_ml.py

#out_nhyd99_mae_complex_4im0_29_1EJ1_18_pp_model_xx01.pdb

grep "^HETATM" ../process_3/out_nhyd99_mae_complex_${v2}_pp_model_xx01.pdb > temp_lig_${v2}.pdb
obabel -ipdb temp_lig_${v2}.pdb -osmi -O temp_smi_${v2}.smi
$openeye tautomers -in temp_smi_${v2}.smi -out tauto_temp_smi_${v2}.smi -maxtoreturn 1 -warts false -stereo None
cat tauto_temp_smi_${v2}.smi >> collected_smiles.smi

