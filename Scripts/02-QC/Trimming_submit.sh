#! /bin/bash

#SBATCH --job-name=trimming             # Job name to show with squeue
#SBATCH --output=trimming_%j.out        # Output file
#SBATCH --ntasks=88                     # Maximum number of cores to use
#SBATCH --time=1-00:00:00               # Time limit to execute the job
#SBATCH --mem-per-cpu=5G                # Required Memory per core
#SBATCH --cpus-per-task=2               # CPUs assigned per task.
#SBATCH --qos=short                     # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   trimming_submit.sh
#   This program executes the 02-Trimming.sh program in parallel to trim the
#   species sequences using Illumina adapters. For this purpose, the Fastp
#   program is utilized.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 11/04/2022
#   Version: 1.3
#
#******************************************************************************

# Paths
path_in=/storage/ncRNA/Projects/sRNA_project/03-Raw_data
path_out=/home/gonsanan/miRNAs_srp_project/Libraries/01-Trimmed_data_20_25nt
path_adapters=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/02-Trimming/01-Trimming_adapters/all_adapters.fa

# List species directories
species_list=$( ls -d1 $path_in/* )

# Iterate files list
for species in $species_list
do
    ## List species projects
    projects_list=$( ls -d1 $species/* )

    ## Iterate species projects
    for project in $projects_list
    do
        ### Execute 01-Prefetch.sh program
        srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive ./Trimming.sh \
            -p $project \
            -o $path_out \
            -a $path_adapters \
            -n 20 \
            -x 25 &
    done
done
wait

exit 0
