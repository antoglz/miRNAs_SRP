#!/bin/bash

#SBATCH --job-name=ident_miRNAs         # Job name to show with squeue
#SBATCH --output=ident_miRNAs_%j.out    # Output file
#SBATCH --ntasks=1                      # Maximum number of cores to use
#SBATCH --time=01-00:00:00              # Time limit to execute the job
#SBATCH --mem=2G                        # Required Memory
#SBATCH --cpus-per-task=20              # CPUs assigned per task.
#SBATCH --qos=short                     # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   miRNAs_annotation_submit.sh
#
#   The program is designed to identify miRNAs from tables of differentially
#   expressed sequences for multiple species using the miRNAS_annotation.sh
#   programm.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 22/12/2023
#   Version: 2.0 
#
#******************************************************************************

# Input paths
path_in=/home/gonsanan/miRNAs_srp_project/Results/06-Diff_exp_analysis
path_mirbase=/storage/ncRNA/Projects/sRNA_project/05-Databases/miRNAs/Sequence_identification/02-Filtered_databases/miRBase
path_PmiREN=/storage/ncRNA/Projects/sRNA_project/05-Databases/miRNAs/Sequence_identification/02-Filtered_databases/PmiREN
path_sRNAanno=/storage/ncRNA/Projects/sRNA_project/05-Databases/miRNAs/Sequence_identification/02-Filtered_databases/sRNAanno
path_ids_table=/storage/ncRNA/Projects/sRNA_project/00-Species_information/species_id.csv
ea_table=/home/gonsanan/miRNAs_srp_project/Results/05-Exploratory_analysis/ea_summary.csv
mww_pvalue=0.05

# Ouput paths
path_out=/home/gonsanan/miRNAs_srp_project/Results/07-Identification_miRNAs

# Threads
num_threads=20

# Execution 
bash miRNAs_annotation.sh \
     --input $path_in \
     --output $path_out \
     --mirbase $path_mirbase \
     --pmiren $path_PmiREN \
     --srnaanno $path_sRNAanno \
     --species-ids $path_ids_table \
     --threads $num_threads \
     --ea-table $ea_table \
     --mww-pvalue $mww_pvalue

exit 0
