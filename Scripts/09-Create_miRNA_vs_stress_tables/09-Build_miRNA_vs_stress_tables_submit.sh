#!/bin/bash

#SBATCH --job-name=pre-abs-table       	# Job name (showed with squeue)
#SBATCH --output=pre_abs_table_%j.out  	# Standard output and error log
#SBATCH --qos=short                	    # The QoS to submit the job.
#SBATCH --nodes=1                   	# Number of nodes requested
#SBATCH --ntasks=1                  	# Required only 1 task
#SBATCH --cpus-per-task=1           	# Required only 1 cpu
#SBATCH --mem=5G                    	# Job memory request.
#SBATCH --time=00:10:00             	# Time limit days-hrs:min:sec.

#******************************************************************************
#  
#   9-Build_presence_absence_table_exe.sh
#
#   This program executes the 9-Build_presence_absence_table.py program.
# 
#   Author: Antonio Gonzalez Sanchez
#   Date: 13/02/2023
#   Version: 1.0 
#
#******************************************************************************

# Modules
source /home/ggomez/joan/ENTER/bin/activate antonio_tfm

# Execution 
python3 ./09-Build_presence_absence_table.py \
    --path-dea /storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/06-DiffExpAnalysis_res/summary.csv \
    --path-annot /storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/08-miRNAs_grouped_by_family \
    --path-dir-out /storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/09-Presence_absence_table

exit 0