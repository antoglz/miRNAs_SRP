#!/bin/bash

#SBATCH --job-name=miRNAs_fam         	# Job name to show with squeue
#SBATCH --output=miRNAs_fam_%j.out    	# Output file
#SBATCH --ntasks=1                 		# Maximum number of cores to use
#SBATCH --time=00-00:30:00          	# Time limit to execute the job
#SBATCH --mem=1G           				# Required Memory
#SBATCH --cpus-per-task=1          		# CPUs assigned per task.
#SBATCH --qos=short                		# QoS: short,medium,long,long-me

#******************************************************************************
#  
#   08-Group_miRNAs_by_family_exe.sh
#
#   This program executes the 08-Group_miRNAs_by_family.r program.
# 
#   Author: Antonio Gonzalez Sanchez
#   Date: 16/02/2023
#   Version: 1.0 
#
#******************************************************************************

# Modules
source /home/ggomez/joan/ENTER/bin/activate antonio_Renv

# Paths
path_in_dea=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/06-DiffExpAnalysis_res
path_in_annot=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/07-Identification_miRNAs_res
path_out=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/08-miRNAs_grouped_by_family

# Other variables
data_type=mature

# Execution 
Rscript 08-Group_miRNAs_by_family.r \
	--dea $path_in_dea \
	--annotation $path_in_annot \
	--type $data_type \
	--output $path_out
    
exit 0
