#!/bin/bash

#SBATCH --job-name=counts_d         # Job name to show with squeue
#SBATCH --output=counts_d_%j.out    # Output file
#SBATCH --ntasks=25                 # Maximum number of cores to use
#SBATCH --time=07-00:00:00          # Time limit to execute the job
#SBATCH --mem-per-cpu=2G            # Required Memory per core
#SBATCH --cpus-per-task=6           # CPUs assigned per task.
#SBATCH --qos=medium                # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   04-Counts_and_sRNADatabase_submit.sh
#
#   This program executes the 04-Counts_and_sRNADatabase.py program in
#   parallel to generate the absolute counts and Reads Per Million
#   (RPM) tables for a specific project using the trimmed and filtered
#   libraries, also calculating the averages of both types of counts in
#   the different replicates of each condition for the sequences analysed.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 20/09/2023
#   Version: 1.1 
#
#******************************************************************************

# Modules
module load anaconda
source activate sRNA_project

# Paths
path_lib=/home/gonsanan/miRNAs_srp_project/Libraries/02-Trimmed_DepthAndRep_filtered/fasta
path_results=/home/gonsanan/miRNAs_srp_project/Results
path_metadata=/storage/ncRNA/Projects/sRNA_project/02-Metadata
path_rnacentral=/storage/ncRNA/Projects/sRNA_project/05-Databases/miRNAs/Sequence_filtering/01-RNAcentral/rnacentral_plants_filtered.fasta

# List all species paths and names
species_list=$( ls -d1 $path_lib/* )

# Iterate specie list
for species in $species_list
do
    ## List species projects
    projects_list=$( ls -d1 $species/* )
    
    ## Iterate species projects
    for project in $projects_list
    do
        ### Execute Counts_and_sRNADatabase.py program
        srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive python3 sRNA_counts.py \
            --path-project $project \
            --path-results $path_results \
            --threads $SLURM_CPUS_PER_TASK \
            --path-metadata $path_metadata \
            --path-rnacentral $path_rnacentral &
    done
done
wait

exit 0
