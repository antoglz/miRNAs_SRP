#!/bin/bash

#SBATCH --job-name=Filter_lib      	# Job name to show with squeue
#SBATCH --output=Filter_lib_%j.out 	# Output file
#SBATCH --ntasks=52             	# Maximum number of cores to use
#SBATCH --time=1-00:00:00         	# Time limit to execute the job
#SBATCH --mem-per-cpu=2G        	# Required Memory per core
#SBATCH --cpus-per-task=2       	# CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#  03-Filter_by_depth_rep_submit.sh
#
#  This program executes the 03-Filter_by_depth_rep.py program in parallel to
#  filter the libraries of the different species at the same time. In addition,
#  it joins the summary files generated for each specie into a single file.
#
#  Author: Antonio Gonzalez Sanchez
#  Date: 25/04/2022
#  Version: 1.1 
#
#******************************************************************************

# Modules
module load anaconda
source activate tsRNA_project

# Paths
path_in=/storage/ncRNA/Projects/tsRNA_project/Libraries/cleanData
path_out=/storage/ncRNA/Projects/tsRNA_project/Libraries/Trimmed_DepthAndRep_filtered
path_out_sum=/storage/ncRNA/Projects/tsRNA_project/Additional_data/DepthAndRep_filtered
path_metadata_dir=/storage/ncRNA/Projects/TFM_AntonioG/Additional_info/Metaanalysis_miRNA/01-Sra_info/Metaanalysis_miRNA/Standardised_metadata

# List all species paths and names
species_list=$( ls -d1 $path_in/* )

# Iterate specie list
for species in $species_list
do
    ## List species projects
    projects_list=$( ls -d1 $species/* )
    
    ## Iterate species projects
    for project in $projects_list
    do
        ### Execute 03-FilterByDepthAndRep.sh program
        srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive python3 Filter_by_depth_rep.py \
                -i $project \
                -o $path_out \
                -s $path_out_sum \
                -m $path_metadata_dir \
                -d 500000 \
                -r 2 \
                -p $SLURM_CPUS_PER_TASK &
    done
wait

# Concatenate files with the results of each specie
awk '{ print $1,$2,$3,$4}' $path_out_sum/*_sum_projects.tsv > $path_out_sum/sum_projects_temp.tsv
sort -t$'\t' -k1,1 -k2  -n $path_out_sum/sum_projects_temp.tsv > $path_out_sum/sum_projects.tsv
sed -i '1s/^/Specie Project Filtered Excluded\n/' $path_out_sum/sum_projects.tsv

# Concatenate summary files of filtered and excluded libraries
[ `ls -1 $path_out_sum/*_sum_libraries.tsv 2>/dev/null | wc -l ` -gt 0 ] && awk '{print}' $path_out_sum/*_sum_libraries.tsv > $path_out_sum/sum_libraries_temp.tsv
sort -t$'\t' -k1 -n $path_out_sum/sum_libraries_temp.tsv > $path_out_sum/sum_libraries.tsv

# Remove specie files
rm -f $path_out_sum/*_sum_results.tsv
rm -f $path_out_sum/*_sum_libraries.tsv
rm -f $path_out_sum/*_temp.tsv

exit 0