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
#SBATCH -J run_sims50
#! Which project should be charged (NB Peta4-KNL projects end in '-KNL'):
#SBATCH -A FUSCO-SL3-KNL
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! Which memory mode do you require (recommended choices are cache or flat):
#SBATCH --constraint=cache
#! How many (MPI) tasks will there be in total? (<=nodes*256)
#SBATCH --ntasks=64
#! How much wallclock time will be required?
#SBATCH --time=02:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! Do not change:
#SBATCH -p knl

#! sbatch directives end here (put any additional directives above this line)

#! Notes:
#! Charging is determined by node number*walltime. Allocation is in entire nodes.

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

sh run_sims_Kb50.sh 
