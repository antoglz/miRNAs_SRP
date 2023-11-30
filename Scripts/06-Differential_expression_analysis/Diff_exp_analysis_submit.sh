#!/bin/bash

#SBATCH --job-name=diffea           # Job name to show with squeue
#SBATCH --output=diffea_%j.out      # Output file
#SBATCH --ntasks=50                 # Maximum number of cores to use
#SBATCH --time=01-00:00:00          # Time limit to execute the job
#SBATCH --mem-per-cpu=2G            # Required Memory per core
#SBATCH --cpus-per-task=2           # CPUs assigned per task.
#SBATCH --qos=short                 # QoS: short,medium,long,long-mem

#******************************************************************************
#  
#   Diff_exp_analysis_submit.sh
#
#   This program executes the Diff_exp_analysis.r program in parallel to
#   perform differential expression analysis using the absolute count tables
#   of each project simultaneously.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 16/02/2021
#   Version: 2.0 
#
#******************************************************************************

# Modules
module load anaconda
source activate sRNA_project

# Paths
path_in=/home/gonsanan/miRNAs_srp_project/Results/04-Projects_divided_by_subprojects
path_ea=/home/gonsanan/miRNAs_srp_project/Results/Metaanalysis_miRNA/05-PCA
path_dea=/home/gonsanan/miRNAs_srp_project/Results/06-DiffExpAnalysis
alpha=0.05

################################################################################
#           1. EXECUTE Diff_exp_analysis.r PROGRAM IN PARALLEL                 #
################################################################################

# List all species paths
species=$( ls -d1 $path_in/* )

# Iterate specie paths list
for path_species in $species
do
    ## List species projects paths
    projects=$( ls -d1 $path_species/* )
    
    ## Iterate species projects paths list
    for path_project in $projects
    do
        ### Execution 
        srun -N1 -n1 -c$SLURM_CPUS_PER_TASK --quiet --exclusive Rscript Diff_exp_analysis.r \
            --input $path_project \
            --output $path_dea \
            --exploratory $path_ea \
            --alpha $alpha &
    done
done
wait

exit 0


################################################################################
#               2. CREATE EXPLORATORY ANALYSIS SUMMARY FILE                    #
################################################################################

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


################################################################################
#           3. CREATE DIFFERENTIAL EXPRESSION ANALYSIS SUMMARY FILE            #
################################################################################

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
