#!/bin/bash

#SBATCH --job-name=pre-abs-table       		# Job name (showed with squeue)
#SBATCH --output=mirnas_vs_stress_table_%j.out  # Standard output and error log
#SBATCH --qos=short                	    	# The QoS to submit the job.
#SBATCH --nodes=1                   		# Number of nodes requested
#SBATCH --ntasks=1                  		# Required only 1 task
#SBATCH --cpus-per-task=1           		# Required only 1 cpu
#SBATCH --mem=5G                    		# Job memory request.
#SBATCH --time=00:10:00             		# Time limit days-hrs:min:sec.

#******************************************************************************
#  
#   Build_miRNA_vs_stress_tables_submit.sh
#
#   This program executes the Build_miRNA_vs_stress_tables.py program.
# 
#   Author: Antonio Gonzalez Sanchez
#   Date: 10/26/2023
#   Version: 1.0 
#
#******************************************************************************

# Modules
module load anaconda
source activate sRNA_project

# Paths
path_dea=/home/gonsanan/miRNAs_srp_project/Results/06-Diff_exp_analysis/summary.csv
path_annot=/home/gonsanan/miRNAs_srp_project/Results/08-miRNAs_grouped_by_family/Group_miRNAs_sig
path_dir_out=/home/gonsanan/miRNAs_srp_project/Results/09-miRNA_vs_stress_tables

# Execution 
python3 ./Build_miRNA_vs_stress_tables.py \
    --path-dea $path_dea \
    --path-annot $path_annot \
    --path-dir-out $path_dir_out

exit 0
