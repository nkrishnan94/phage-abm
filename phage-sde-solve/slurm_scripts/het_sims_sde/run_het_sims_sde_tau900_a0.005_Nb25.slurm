#!/bin/bash
#SBATCH -J run_sims_Nb100tau900_a0.005_Nb25
#SBATCH -A FUSCO-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=5980mb
#SBATCH --array=1-500
#SBATCH -p skylake
./etc/profile.d/modules.sh
module purge 
module load rhel7/default-peta4
echo "This is job" $SLURM_ARRAY_TASK_ID
g++ -o outanc_tau900_a0.005_Nb25_$SLURM_ARRAY_TASK_ID phage_inf_coarse_het_Nb25.cpp -std=c++11
./outanc_tau900_a0.005_Nb25_$SLURM_ARRAY_TASK_ID -t 900 -a 0.005 -i $SLURM_ARRAY_TASK_ID
