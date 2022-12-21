#!/bin/bash
#SBATCH --cpus-per-task=60
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=0               # memory per node
#SBATCH --job-name=omega

cp=60
$openeye oeomega classic -in $1 -out $1.sdf  -flipper $2 -strictstereo false -maxconfs 1 -mpi_np $cp -log $1.log -prefix $1
