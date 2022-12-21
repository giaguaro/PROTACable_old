#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=results

source ~/.bashrc
conda activate py39

cd $1
original=$PWD

cd process_4

rm ${original}/ml_scores/*pdb
rm ${original}/ml_scores/*all_ml_scores.csv
rm ${original}/ml_scores/*smi
rm ${original}/ml_scores/collected_smiles.smi


for p in *xx01-out-out-1.pdb; do sbatch ml_scores.sh $p;done

while [[ $(squeue --name scores | wc -l) -gt 1 ]]; do echo "wait"; sleep 1s; done

for i in {1..9}; do for p in *xx0${i}-out-out-1.pdb; do id=$(echo "$p" | cut -d '_' -f6-9); if [ ! -f ../ml_scores/${id}_all_ml_scores.csv ]; then sbatch --wait ml_scores.sh $p;fi;done;done

for i in {10..20}; do for p in *xx${i}-out-out-1.pdb; do id=$(echo "$p" | cut -d '_' -f6-9); if [ ! -f ../ml_scores/${id}_all_ml_scores.csv ]; then sbatch --wait ml_scores.sh $p;fi;done;done

cd ${original}/ml_scores


python ml_learning_dummy.py

mkdir results
rm results/*

while read line; do for r in ../process_4/*${line}*csv; do f1=$(awk -F',' 'NR==3 {print $2}' $r); v1=${r%%-*}; vf=$(echo ${v1} | cut -d '/' -f 3| cut -d '_' -f 2-12); echo "${vf} ${f1}" >> temp_chosen.txt; done; c=$(sort -k2 -n temp_chosen.txt| head -n1 | cut -d ' ' -f1); echo $c; rm temp_chosen.txt; cp ../process_3/${c}.pdb results; done <chosen.txt

for i in results/*; do sed -i '/HETATM/s/  BR  / BR   /g' $i;done
for i in results/*; do sed -i '/HETATM/s/  CL  / CL   /g' $i;done

