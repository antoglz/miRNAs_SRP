#!/bin/bash

#SBATCH --job-name=prefetch_lib         # Job name to show with squeue
#SBATCH --output=prefetch_lib_%j.out    # Output file
#SBATCH --ntasks=48                     # Maximum number of cores to use
#SBATCH --time=1-00:00:00               # Time limit to execute the job
#SBATCH --mem-per-cpu=5G                # Required Memory per core
#SBATCH --cpus-per-task=1               # CPUs assigned per task.
#SBATCH --qos=short                     # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   Prefetch_submit.sh
#
#   This program executes the Prefetch.sh program in parallel to
#   download the libraries of a group of SRA projects and perform a
#   quality control check (FastQC).
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 05/04/2022
#   Version: 1.0 
#
#******************************************************************************

# Paths
path_in=/storage/ncRNA/Projects/tsRNA_miRNA_project/01-Accession_list
path_out=/storage/ncRNA/Projects/tsRNA_miRNA_project/03-Raw_data

# List project files with accession list
project_acc_path_list=$( ls -d1 $path_in/* )

# Iterate files list
for project_acc_path in $project_acc_path_list
do
    ## Execute Prefetch.sh program
    srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive ./Prefetch.sh \
        -p $project_acc_path \
        -o $path_out &
done
wait

exit 0
