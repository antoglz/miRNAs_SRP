#!/bin/bash

#SBATCH --job-name=miRNAs_fam         		# Job name to show with squeue
#SBATCH --output=miRNAs_fam_%j.out    		# Output file
#SBATCH --ntasks=1                 		# Maximum number of cores to use
#SBATCH --time=00-00:30:00          		# Time limit to execute the job
#SBATCH --mem=1G           			# Required Memory
#SBATCH --cpus-per-task=1          		# CPUs assigned per task.
#SBATCH --qos=short                		# QoS: short,medium,long,long-me

#******************************************************************************
#  
#   Group_miRNAs_by_family_submit.sh
#
#   This program executes the Group_miRNAs_by_family.r program.
# 
#   Author: Antonio Gonzalez Sanchez
#   Date: 12/20/2023
#   Version: 2.0 
#
#******************************************************************************

# Modules
module load anaconda
source activate sRNA_project

# Paths
path_in_dea=/home/gonsanan/miRNAs_srp_project/Results/06-Diff_exp_analysis
path_in_annot=/home/gonsanan/miRNAs_srp_project/Results/07-Identification_miRNAs
path_out=/home/gonsanan/miRNAs_srp_project/Results/08-miRNAs_grouped_by_family

# Other variables
data_type=mature

# Execution 
Rscript Group_miRNAs_by_family.r \
	--dea $path_in_dea \
	--annotation $path_in_annot \
	--type $data_type \
	--output $path_out
    
exit 0
