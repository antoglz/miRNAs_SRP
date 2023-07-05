#! /bin/bash

#SBATCH --job-name=prefetch_lib         # Job name to show with squeue
#SBATCH --output=prefetch_lib_%j.out    # Output file
#SBATCH --ntasks=48                     # Maximum number of cores to use
#SBATCH --time=1-00:00:00               # Time limit to execute the job
#SBATCH --mem-per-cpu=5G                # Required Memory per core
#SBATCH --cpus-per-task=1               # CPUs assigned per task.
#SBATCH --qos=short                     # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   01-Prefetch_parallel.sh
#
#   This program executes the 01-Prefetch.sh program in parallel to
#   download the libraries of a group of SRA projects and perform a
#   quality control check.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 05/04/2022
#   Version: 1.0 
#
#******************************************************************************

## Paths
path_in=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/01-Sra_info/Metaanalysis_miRNA/Accession_list
path_out=/storage/ncRNA/Projects/TFM_AntonioG/Libraries/Metaanalysis_miRNA/01-Raw_data/Illumina

## List project files with accession list
project_path_list=$( ls -d1 $path_in/* )

## Iterate files list
for project_path in $project_path_list
do
    # Execute 01-Prefetch.sh program
    srun -N1 -n1 -c1 --quiet --exclusive ./01-Prefetch.sh \
        -p $project_path \
        -o $path_out &
done
wait

exit 0
