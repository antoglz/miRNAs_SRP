#!/bin/bash

#SBATCH --job-name=divide_prj           # Job name to show with squeue
#SBATCH --output=divide_prj_%j.out      # Output file
#SBATCH --ntasks=30                     # Maximum number of cores to use
#SBATCH --time=01-00:00:00              # Time limit to execute the job
#SBATCH --mem-per-cpu=7G                # Required Memory per core
#SBATCH --cpus-per-task=2               # CPUs assigned per task.
#SBATCH --qos=short                     # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   05-Divide_filter_by_subproject_submit.sh
#
#   This program executes the 05-Divide_filter_by_subprojects.py
#   program in parallel to simultaneously divide by experiment the count
#   tables of the different projects, writing the table section of each
#   experiment in a separate file.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 21/09/2023
#   Version: 1.1 
#
#******************************************************************************

# Modules
module load anaconda
source activate tsRNA_project

# Paths
path_counts=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/03-Fusion_count_tables_RF
metadata=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/01-Sra_info/Metaanalysis_miRNA/Standardised_metadata
path_results=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/04-Projects_divided_by_subprojects

# List all species paths
species_list=$( ls -d1 $path_counts/* )

# Iterate specie paths list
for species in $species_list
do
    ## List species projects paths
    projects_list=$( ls -d1 $species/* )
    
    ## Iterate species projects paths list
    for project in $projects_list
    do
        ### Execution
        srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive python3 ./Divide_filter_by_subprojects.py \
            --path-counts $sp \
            --metadata $metadata \
            --path-results $path_results &
    done
done
wait

exit 0

