#!/bin/bash

#SBATCH --job-name=ident_miRNAs         # Job name to show with squeue
#SBATCH --output=ident_miRNAs_%j.out    # Output file
#SBATCH --ntasks=1                 	    # Maximum number of cores to use
#SBATCH --time=01-00:00:00          	# Time limit to execute the job
#SBATCH --mem=2G           		        # Required Memory
#SBATCH --cpus-per-task=30          	# CPUs assigned per task.
#SBATCH --qos=short                	    # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   07-miRNAs_identification_exe.sh
#
#   The program is designed to identify miRNAs from tables of differentially
#   expressed sequences for multiple species using the 07-miRNAS_identifi-
#   cation.sh programm.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 16/02/2021
#   Version: 1.0 
#
#******************************************************************************

# Input paths
path_in=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/06-DiffExpAnalysis_res
path_mirbase=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/10-miRNAsAnnotation/01-Databases/miRBase
path_PmiREN=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/10-miRNAsAnnotation/01-Databases/PmiREN
path_ids_table=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/10-miRNAsAnnotation/species_id.csv
mww_pvalue=0.05
ea_table=/home/antonio/Escritorio/pruebas/prueba_tabla_anotación/Results/05-PCA/ea_summary.csv

# Ouput paths
path_out=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/07-Identification_miRNAs_res

# Threads
num_threads=30

# Execution 
bash miRNAs_annotation.sh \
    --input $path_in \
    --output $path_out \
    --mirbase $path_mirbase \
    --pmiren $path_PmiREN \
    --species-ids $path_ids_table \
    --threads $num_threads \
    --mww-pvalue $mww_pvalue \
    --ea-table $ea_table


exit 0
