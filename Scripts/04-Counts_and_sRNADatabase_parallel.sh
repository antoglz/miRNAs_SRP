#!/bin/bash

#SBATCH --job-name=counts_d         # Job name to show with squeue
#SBATCH --output=counts_d_%j.out    # Output file
#SBATCH --ntasks=15                 # Maximum number of cores to use
#SBATCH --time=07-00:00:00          # Time limit to execute the job
#SBATCH --mem-per-cpu=2G            # Required Memory per core
#SBATCH --cpus-per-task=10          # CPUs assigned per task.
#SBATCH --qos=medium                # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   04-Counts_and_sRNADatabase_parallel.sh
#
#   This program executes the 04-Counts_and_sRNADatabase.py program in
#   parallel to generate the absolute counts and Reads Per Million (RPMs)
#   tables fo the different species using the trimmed and filtered libraries,
#   also calculating the averages of both types of counts in the different
#   replicates of each condition for the sequences analysed. 
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 16/02/2021
#   Version: 1.0 
#
#******************************************************************************

# Modules
source /home/ggomez/joan/ENTER/bin/activate antonio_tfm

# Paths
path_lib=/storage/ncRNA/Projects/TFM_AntonioG/Libraries/Metaanalysis_miRNA/03-Trimmed_DepthAndRep_filtered
path_results=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA
path_database=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/05-sRNADatabases
path_metadata=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/01-Sra_info/Metaanalysis_miRNA/Standardised_metadata
path_rnacentral=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/04-RNAcentral/rnacentral_plants_filtered.fasta

# Other variables
create_database=yes
max_digits=10
threads=10

# List all species paths and names
species=$( ls -d1 $path_lib/* )

# Iterate species list
for sp in $species
do
    ## Execution 
    srun -N1 -n1 -c10 --quiet --exclusive python3 04-Counts_and_sRNADatabase.py \
        --path-species $sp \
        --path-results $path_results \
        --create-database $create_database \
        --path-database $path_database \
        --max-digits $max_digits \
        --threads $threads \
        --path-metadata $path_metadata \
        --path-rnacentral $path_rnacentral &
done
wait

exit 0
