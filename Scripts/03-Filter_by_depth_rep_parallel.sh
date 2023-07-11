#!/bin/bash

#SBATCH --job-name=Filter_lib      	# Job name to show with squeue
#SBATCH --output=Filter_lib_%j.out 	# Output file
#SBATCH --ntasks=32             	# Maximum number of cores to use
#SBATCH --time=15:00:00         	# Time limit to execute the job
#SBATCH --mem-per-cpu=1G        	# Required Memory per core
#SBATCH --cpus-per-task=1       	# CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#  03-Filter_by_depth_rep_parallel.sh
#
#  This program executes the 03-Filter_by_depth_rep.py program in parallel to
#  filter the libraries of the different species at the same time. In addition,
#  it joins the summary files generated for each specie into a single file.
#
#  Author: Antonio Gonzalez Sanchez
#  Date: 25/04/2022
#  Version: 1.0 
#
#******************************************************************************

# Modules
source /home/ggomez/joan/ENTER/bin/activate antonio_tfm

# Paths
path_in="/storage/ncRNA/Projects/TFM_AntonioG/Libraries/Metaanalysis_miRNA/02-Trimmed_data"
path_out_fasta="/storage/ncRNA/Projects/TFM_AntonioG/Libraries/Metaanalysis_miRNA/03-Trimmed_DepthAndRep_filtered"
path_out_sum="/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/03-Filtered_projects/03-DepthAndRep_filtered"
path_metadata_dir="/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/01-Sra_info/Metaanalysis_miRNA/Standardised_metadata"

# List all species paths and names
species=$( ls -d1 $path_in/* )

# Iterate specie list
for sp in $species
do
    # Execute 03-FilterByDepthAndRep.sh program
    srun -N1 -n1 -c1 --quiet --exclusive python3 03-Filter_by_depth_rep.py \
        -i $sp \
        -o $path_out_fasta \
        -s $path_out_sum \
        -m $path_metadata_dir \
        -d 500000 \
        -r 2 &
done
wait

# Concatenate files with the results of each specie
awk '{ print $1,$2,$3,$4}' $path_out_sum/*_results.txt > $path_out_sum/results.txt
sed -i '1s/^/Specie Project Filtered Excluded\n/' $path_out_sum/results.txt

# Concatenate summary files of filtered and excluded libraries
[ `ls -1 $path_out_sum/*_filtered.txt 2>/dev/null | wc -l ` -gt 0 ] && awk '{print}' $path_out_sum/*_filtered.txt > $path_out_sum/filtered.txt
[ `ls -1 $path_out_sum/*_excluded.txt 2>/dev/null | wc -l ` -gt 0 ] && awk '{print}' $path_out_sum/*_excluded.txt > $path_out_sum/excluded.txt

# Remove specie files
rm -f $path_out_sum/*_results.txt
rm -f $path_out_sum/*_filtered.txt
rm -f $path_out_sum/*_excluded.txt

exit 0


