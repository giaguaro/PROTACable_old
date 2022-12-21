#!/bin/bash
#SBATCH --time=168-00:00
#SBATCH --partition=normal
#SBATCH --job-name=icm
#SBATCH --ntasks=1
#SBATCH --array=0-1                 #IMPORTANT: do not go over 1000 for y in x-y otherwise it will not run. Modify size_single below instead. y*$size_single must be >= tot_compounds-$size_single-1

###VARIABLES
project=$1
input_sdf=$2
thoroughness=1.0
output=$3
size_single=2                              #compounds in each job (first job will have $size_single-1)
in=$((SLURM_ARRAY_TASK_ID*size_single))
end=$((in+size_single-1))

###ICM RUN
$ICMHOME/icm -vlscluster $ICMHOME/_dockScan $project -a input=$input_sdf -S confs=10 maxFileSizeMb=10000 effort=$thoroughness from=$in to=$end name=$output'_'$in output=$output'_'$in >> $output'_'$in'.'log
