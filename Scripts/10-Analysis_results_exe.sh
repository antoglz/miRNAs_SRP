#!/bin/bash

#SBATCH --job-name=results           # Job name to show with squeue
#SBATCH --output=results_%j.out      # Output file
#SBATCH --ntasks=1                 # Maximum number of cores to use
#SBATCH --time=00-1:00:00          # Time limit to execute the job
#SBATCH --mem-per-cpu=2G            # Required Memory per core
#SBATCH --cpus-per-task=1           # CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   10-Analysis_results_exe.sh
#
#   This program executes the 10-Analysis_results.r program.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 16/02/2021
#   Version: 1.0 
#
#******************************************************************************

# Modules
source /home/ggomez/joan/ENTER/bin/activate antonio_Renv

# Paths
pre_abs_table_path=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/09-Presence_absence_table/presence_absence_table_sig.csv
ids_path=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/09-Presence_absence_table/ids_table.csv
species_name_path=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/10-miRNAsAnnotation/species_id.csv
sRNAs_sig_path=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/06-DiffExpAnalysis_res/summary.csv
miRNAs_ident_path=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/07-Identification_miRNAs_res/miRNA_identification_sig/summary.csv
output_directory_path=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/10-Analysis_results

# Execute 10-Analysis_results.r program
Rscript 10-Analysis_results.r \
    --preabstable $pre_abs_table_path \
    --ids $ids_path \
    --spnames $species_name_path \
    --diffexp $sRNAs_sig_path \
    --mirnasannot $miRNAs_ident_path \
    --output $output_directory_path

exit 0
