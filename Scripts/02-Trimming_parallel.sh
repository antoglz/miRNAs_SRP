#! /bin/bash

#SBATCH --job-name=trimming         # Job name to show with squeue
#SBATCH --output=trimming_%j.out    # Output file
#SBATCH --ntasks=32                 # Maximum number of cores to use
#SBATCH --time=1-00:00:00           # Time limit to execute the job
#SBATCH --mem-per-cpu=5G            # Required Memory per core
#SBATCH --cpus-per-task=1           # CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   02-Trimming_parallel.sh
#
#   This program executes the 02-Trimming.sh program in parallel to trim the
#   species sequences using Illumina adapters. For this purpose, the Fastp
#   program is utilized.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 11/04/2022
#   Version: 1.0 
#
#******************************************************************************

# Paths
path_in=/storage/ncRNA/Projects/TFM_AntonioG/Libraries/Metaanalysis_miRNA/01-Raw_data/Illumina
path_out=/storage/ncRNA/Projects/TFM_AntonioG/Libraries/Metaanalysis_miRNA/02-Trimmed_data
path_adapters=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/02-Trimming/01-Trimming_adapters/all_adapters.fa

# List speciess directories
species_list=$( ls -d1 $path_in/* )

# Iterate files list
for species in $species_list
do
    ## Execute 01-Prefetch.sh program
    srun -N1 -n1 -c1 --quiet --exclusive ./02-Trimming.sh \
        -s $species \
        -o $path_out \
        -a $path_adapters &
done
wait

exit 0
