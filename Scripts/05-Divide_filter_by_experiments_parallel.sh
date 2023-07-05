#!/bin/bash

#SBATCH --job-name=divide_prj           # Job name to show with squeue
#SBATCH --output=divide_prj_%j.out      # Output file
#SBATCH --ntasks=30                     # Maximum number of cores to use
#SBATCH --time=01-00:00:00              # Time limit to execute the job
#SBATCH --mem-per-cpu=7G                # Required Memory per core
#SBATCH --cpus-per-task=1               # CPUs assigned per task.
#SBATCH --qos=short                     # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   05-Divide_filter_by_experiments_parallel.sh
#
#   This program executes the 05-Divide_filter_by_experiments.py
#   program in parallel to simultaneously divide by experiment the count
#   tables of the different species, writing the table section of each
#   experiment in a separate file.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 21/12/2022
#   Version: 1.0 
#
#******************************************************************************

# Modules
source /home/ggomez/joan/ENTER/bin/activate antonio_tfm

# Paths
path_counts=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/03-Fusion_count_tables_RF
metadata=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/01-Sra_info/Metaanalysis_miRNA/Standardised_metadata
path_results=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/04-Projects_divided_by_experiments

# List all species paths and names
species=$( ls -d1 $path_counts/* )

# Iterate specie list
for sp in $species
do
    ## Execution
    srun -N1 -n1 -c1 --quiet --exclusive python3 ./05-Divide_filter_by_experiments.py \
        --path-counts $sp \
        --metadata $metadata \
        --path-results $path_results &
done
wait

exit 0
