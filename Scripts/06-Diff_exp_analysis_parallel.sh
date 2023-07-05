#!/bin/bash

#SBATCH --job-name=diffea           # Job name to show with squeue
#SBATCH --output=diffea_%j.out      # Output file
#SBATCH --ntasks=29                 # Maximum number of cores to use
#SBATCH --time=00-20:00:00          # Time limit to execute the job
#SBATCH --mem-per-cpu=2G            # Required Memory per core
#SBATCH --cpus-per-task=1           # CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   06-Diff_exp_analysis_parallel.sh
#
#   This program executes the 06-DiffExpAnalysis.r program in parallel to
#   perform differential expression analysis using the absolute count tables
#   of each species simultaneously.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 16/02/2021
#   Version: 1.0 
#
#******************************************************************************

# Modules
source /home/ggomez/joan/ENTER/bin/activate antonio_Renv

# Paths
path_in=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/04-Projects_divided_by_experiments
path_dea=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/06-DiffExpAnalysis_res
path_ea=/storage/ncRNA/Projects/TFM_AntonioG/Results/Metaanalysis_miRNA/05-PCA
alpha=0.05

#### 1. Execute 09-DiffExpAnalysis.r program in parallel

# List all species paths and names
species=$( ls -d1 $path_in/* )

# Iterate species list
for path_species in $species
do
    ## Execution 
    srun -N1 -n1 -c1 --quiet --exclusive Rscript 06-Diff_exp_analysis.r \
        --input $path_species \
        --output $path_dea \
        --exploratory $path_ea \
        --alpha $alpha &
done
wait


#### 2. Create exploratory analysis summary file

# If there is more than one file...
if [[ $(ls $path_ea/*_ea_summary_table.csv | wc -l) > 1 ]]
then
    # Concatenate exporatory analysis summary files
    cat $(ls $path_ea/*_ea_summary_table.csv | head -n1) > $path_ea/ea_summary.csv && tail -n +2 -q $path_ea/*_ea_summary_table.csv >> $path_ea/ea_summary.csv
    # Delete individual species files
    rm -f $path_ea/*_ea_summary_table.csv
# If there is only one file...
else
    # Rename it
    mv $path_ea/*_ea_summary_table.csv $path_ea/ea_summary.csv
fi


#### 3. Create differential expression analysis summary file

# If there is more than one file...
if [[ $(ls $path_dea/*_summary.csv | wc -l) > 1 ]]
then
    # Concatenate differential expression analysis summary files
    cat $(ls $path_dea/*_summary.csv | head -n1) > $path_dea/summary.csv && tail -n +2 -q $(ls $path_dea/*_summary.csv | tail -n+2) >> $path_dea/summary.csv
    # Delete individual species files
    rm -f $path_dea/*_summary.csv
# If there is only one file...
else
    # Rename it
    mv $path_dea/*_summary.csv $path_dea/summary.csv
fi

exit 0
