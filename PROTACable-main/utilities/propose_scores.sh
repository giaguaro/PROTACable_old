#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=PROTACs

source ~/.bashrc
conda activate py39


rm *scores pp_scores.csv
awk '/SCORES:/,/Saved/' $1  >> trial.scores
tail -n +2 trial.scores >> trial2.scores
head -n -1 trial2.scores > temp.txt ; mv temp.txt trial2.scores
column -t trial2.scores >> trial3.scores
sed -e 's/\s\+/,/g' trial3.scores > pp_scores.csv
#sed -i '2,102s/^.//' pp_scores.csv
rm *scores
