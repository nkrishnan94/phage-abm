#!/bin/bash
#!
#! Example SLURM job script for Peta4-KNL (KNL 7210)
#! Last updated: Mon 13 Nov 12:17:14 GMT 2017
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J run_sims_test
#! Which project should be charged (NB Peta4-KNL projects end in '-KNL'):
#SBATCH -A FUSCO-SL3-KNL
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! Which memory mode do you require (recommended choices are cache or flat):
#SBATCH --constraint=cache
#! How many (MPI) tasks will there be in total? (<=nodes*256)
#SBATCH --ntasks=64
#! How much wallclock time will be required?
#SBATCH --time=12:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#! Submit a job array with index values between 0 and 31
#! NOTE: This must be a range, not a single number (i.e. specifying '32' here would only run one job, with index 32)
#SBATCH --array=1-100

#! This is the partition name.
#SBATCH -p knl

#! mail alert at start, end and abortion of execution
#! emails will default to going to your email address
#! you can specify a different email address manually if needed.
##SBATCH --mail-type=ALL

#! Don't put any #SBATCH directives below this line

#! Modify the environment seen by the application. For this example we need the default modules.
. /etc/profile.d/modules.sh                # This line enables the module command
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#! The variable $SLURM_ARRAY_TASK_ID contains the array index for each job.
#! In this example, each job will be passed its index, so each output file will contain a different value
echo "This is job" $SLURM_ARRAY_TASK_ID

#! Command line that we want to run:
g++ -o outfp5_$SLURM_ARRAY_TASK_ID phage_inf_sim_het.cpp -std=c++11 -lgsl 
./outfp5_$SLURM_ARRAY_TASK_ID -t .5 -i $SLURM_ARRAY_TASK_ID
