#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#******************************************************************************
#  
#   Filter_by_depth_rep.py
#
#   This program filters the libraries of a given project by selecting those 
#   that exceed the minimum thresholds of sequencing depth and number of
#   replicates specified by the user. The selected libraries will be stored
#   in FASTA and FASTQ format in a output directory specified by the user. 
#   Additionally,this script generates 3 files: species_filtered.txt, which
#   contains the absolute paths of the libraries that exceed the thresholds
#   of sequencing depth and number of replicates; species_excluded.txt,
#   which contains the absolute paths of the libraries that do not exceed
#   the previously mentioned thresholds; and species_results.txt, where the
#   number of libraries that exceed and do not exceed the thresholds in each
#   project is represented. This is represented as a table, where the first
#   column corresponds to the species, the second one to the project, the
#   third one to the number of libraries that exceed the thresholds and the
#   fourth one to the number of libraries that do not exceed it.
#
#   Author: Antonio Gonzalez Sanchez and Julia Corell
#   Date: 26/07/2022
#   Version: 2.0 
#
#
#******************************************************************************


### IMPORT MODULES

import argparse
from Bio import SeqIO
import numpy as np
import multiprocessing
import os
import pandas as pd
import subprocess
import sys


### FUNCTIONS

def GetProjectMetadata(path: str) -> pd.DataFrame:
    '''
    This function extracts the sample metadata from the metadata table of the
    project to which it belongs. 

    Parameters
    ----------
    path : str
        Absolute path of the txt file containing the project metadata.

    Returns
    -------
    df : pd.DataFrame
        Metadata dataframe of the sample of interest.
    '''

    try:
        # Read metadata csv file
        metadata_table = np.genfromtxt(path, delimiter=',', dtype=None, encoding=None,  names=True)

        # Create dataframe
        df = pd.DataFrame(metadata_table)

    except Exception as e:
        print('Unable to read metadata table.')
        print('Exception: ', e)
        sys.exit()
    else:
        return df


def GetInformativeVariables(df: pd.DataFrame) -> list:
    '''   
    This function identifies informative variables within a metadata table,
    understanding informative variables as those with more than one level
    (e.g. Time: 0h, 3h and 6h).
    
    Parameters
    ----------
    df : pd.Dataframe
        Metadata table

    Returns
    -------
    inf_var_pos : list
        List of column indices corresponding to informative variables within
        a list of metadata dataframe columns
    '''
    # List of columns or indexes in which the informative variables are
    # found within the data.frame. Indexes 3, 4 and 5 are condition
    # (cntl vs trtd), stress (p.e. cold) and replicate (p.e. 1,2 or 3),
    # which will always be considered.
    inf_var_pos = [3, 4, 5]
    
    # Get the number of columns except the first 6 (already saved or
    # not important).
    ncol = len(df.columns) - 6
    
    # Iterate ncol
    for col in range(ncol):
        # Extract the "col" column
        #col_list =  list(df.iloc[:,col + 6])
        col_series =  df.iloc[1:,col + 6]

        # Check if the variable contains different levels
        #niveles = pd.unique(col_list[1:])
        levels = col_series.unique()

        # Save the indexes of those variables that have more than 1 level.
        if len(levels) != 1 :
            indice = col + 6
            inf_var_pos.append(indice)

    return inf_var_pos


def GetLibrariesDepthPROCCESSING(path_project: str, path_project_metadata: str,
                                 depth_threshold: int, num_process: int):
    '''
    This function is used to parallelize the GetLibrariesDept function, which
    obtains the sequencing depth of each library within a group of libraries.

    Parameters
    ----------
    path_library_list : list
        List with the paths of compressed FASTQ files (.gz) from a group
        of libraries.
    path_project_metadata : str
        Path to the project's metadata table to which the libraries to be
        used belong.
    depth_threshold: int
        Sequencing depth threshold.
    num_process : int
        Number of processes.

    '''

    # Get a list of libraries and their paths
    libraries_list = os.listdir(path_project)
    libraries_path_list = [f'{path_project}/{library}' for library in libraries_list]

    # Create shared dictionary
    manager = multiprocessing.Manager()
    depth_info_lib = manager.dict() 


    # Calculate process distribution
    number_files = len(libraries_list)
    # If the number of libraries is <= to the max number of processes...
    if number_files <= num_process:
        # Number of processes = number of libraries
        num_process = number_files
        # And, therefore, 1 library for each process
        distribution = 1
    # If there are more than 30 libraries, num_process will always be the
    # maximum (30 in this case).
    else:
        # The libraries in the list are distributed among 30 processes
        distribution = number_files // num_process

    # "pi" would be the starting point of the list section assigned to each
    # process and "pf" the end point.
    processes = [] 
    pi = 0
    for i in range(num_process):
        
        if i == num_process - 1:
            pf = number_files
        else:
            pf = pi + distribution
        
        processes.append(multiprocessing.Process(target = GetLibrariesDepth,
                                                 args=(libraries_path_list[pi:pf],
                                                       path_project_metadata,
                                                       depth_threshold,
                                                       depth_info_lib) )) 
        processes[i].start()
        print(f'Process {i} launched.')
        
        pi = pf
        
    for p in processes:
        p.join()
    
    return(depth_info_lib)


def GetLibrariesDepth(path_library_list: list, path_project_metadata: str,
                      depth_threshold: int, depth_info_lib: multiprocessing.Manager()):
    """
    This function obtains the sequencing depth of each library within a group of
    libraries (understood as the total number of sequences in the library). Once
    obtained, it enters the value into a dictionary shared between processes
    along with a string that indicates whether the library is valid (valid) or
    not (not-valid) based on a depth threshold.

    Parameters
    ----------
    path_library_list : list
        List with the paths of compressed FASTQ files (.gz) from a group
        of libraries.
    path_project_metadata : str
        Path to the project's metadata table to which the libraries to be
        used belong.
    depth_threshold: int
        Sequencing depth threshold.
    depth_info_lib : multiprocessing.Manager()
        Shared dictionary among processes.
    """

    # Iterate libraries
    for path_library in path_library_list:

        # Get library name
        fastq_gz_name = os.path.basename(path_library)

        # Check if the library is in the metadata
        run_name = fastq_gz_name.replace('_tr.fastq.gz','')
        lib_in_metadata = path_project_metadata.loc[path_project_metadata['Run'].str.contains(run_name, case=False)].any().any()
        
        # The library is in the metadata
        if lib_in_metadata:
        
            # Decompress fastq.gz file
            os.system('unpigz -f ' + path_library)

            # Create fastq and fasta file paths
            fastq_path =  path_library.replace('.fastq.gz', '.fastq')
            fasta_path =  fastq_path.replace('.fastq', '.fasta')

            # Convert fastq to fasta
            fastq_record = SeqIO.parse(fastq_path, 'fastq')
            SeqIO.write(fastq_record, fasta_path, 'fasta')

            # Compress fastq file
            os.system('pigz --fast -f ' + fastq_path)

            # Calculate sequencing depth
            cmd = f'grep "^>" {fasta_path} | wc -l'
            depth = int(subprocess.check_output(cmd, shell=True).decode('utf-8').strip())

            # Filter by sequencing depth
            if depth >= depth_threshold:
                depth_info_lib[fastq_gz_name] = [depth, 'valid']
            else:
                depth_info_lib[fastq_gz_name] = [depth, 'not-valid']
        else:
            # Write not-in-metadata) if the sample is not found in the metadata
            depth_info_lib[fastq_gz_name] = ['NA', 'not-in-metadata']
    
    return(depth_info_lib)

    

def GetProjectGroups(list_srr, path_metadata):
    '''
    This function receives a list of samples (SRA Run) from a project and groups
    them according to the experiment they belong to. It does this by comparing
    control samples with treatment samples to identify which treatment samples
    are associated with which control samples. Once identified, it groups them
    in a list together with their control samples. The function returns a list
    containing the lists corresponding to each of these groups, where each
    sample is presented with its Run and its condition separated by "/"
    (e.g. SRR1848795/treated_drought_30%_T.0_leaves_27d).

    Parameters
    ----------
    list_srr : list
        Sample list (SRA Run)
    path_metadata : str
        Metadata table path

    Returns
    -------
    list
        List of control-treated groups
    '''

    # Read metadata and find informative variables
    metadata = GetProjectMetadata(path_metadata)
    inf_var = GetInformativeVariables(metadata)
    inf_var.append(2)

    # Get the name of the columns corresponding to the informative variables.
    inf_colnames = metadata.columns[inf_var].to_list()
    inf_colnames.remove('Replicate')
    inf_colnames.remove('Run')
    
    # Join informative variables in a column (conditions)
    metadata['Condition_all'] = metadata[metadata.columns[inf_var]].apply(lambda x: '_'.join([str(i).strip() for i in x]), axis=1)

    # Extract conditions from dataframe
    samples_list = [condition.replace(' ','').replace(',','_') for condition in metadata['Condition_all'].to_list()]
    
    # Extract metadata from the libraries stored in list_srr
    samples_srr_list = []
    for srr in list_srr:
        for sample_met in samples_list:
            if srr in sample_met:
                samples_srr_list.append(sample_met)

    # Samples list (with replicate number)
    control_samples = []
    treated_samples = []

    # Conditions list (without replicate number)
    control_conditions = []
    treated_conditions = []

    # Separate controls and treatments into two separate lists
    for sample in samples_srr_list:
        # Delete replicate and run columns
        sample_elements = sample.split('_')
        del sample_elements[-1]
        del sample_elements[2]
        # Save control samples
        if sample.startswith('control'):
            # Save control samples in list 
            control_samples.append(sample)
            # Save control conditions in list
            if sample_elements not in control_conditions:
                control_conditions.append(sample_elements)
                
        # Save treated conditions in list
        elif sample.startswith('treated'):
            treated_samples.append(sample)
            if sample_elements not in treated_conditions:
                treated_conditions.append(sample_elements)
    
    # Create matching matrix
    matrix = np.zeros((len(control_conditions), len(treated_conditions)))

    # Check the similarity between controls and treatments
    for i, control_sample in enumerate(control_conditions):
        for j, treated_sample in enumerate(treated_conditions):
            for k, column in enumerate(inf_colnames):
                # If the level is not defined, continue 
                if '.0' in control_sample[k]:
                    continue
                # If there is more than one level, separate them.
                if control_sample[k].find(':'):
                    controrl_sample_ele = control_sample[k].split(':')
                    # If the levels match, add one to the counter.
                    if treated_sample[k] in controrl_sample_ele:
                        matrix[i][j] += 1
                # If the levels match, add one to the counter.
                elif treated_sample[k] == control_sample[k]:
                        matrix[i][j] += 1
    
    # Grouping controls and treatments by experiment
    groups_list = [[] for x in range(len(control_conditions))]
    done = []
    # Iterate matrix rows
    for i in range(len(matrix)):
        # Get maximum value of the row
        max_v = max(matrix[i])
        # Add control samples to the groups_list
        for control in control_samples:
            # Check if it has been previously saved
            if control in done:
                continue
            # Discard replicate and run data
            control_elements = control.split('_')
            run_value = control_elements[-1]
            del control_elements[-1]
            del control_elements[2]
            # Check if conditions match
            if control_elements == control_conditions[i]:
                # Save in group list
                done.append(control)
                con = '_'.join(control_elements)
                groups_list[i].append(f'{run_value}/{con}')
        
        # Add treated samples
        for j in range(len(matrix[i])):
            if matrix[i][j] == max_v:
                # Add treated samples to the groups_list
                for treated in treated_samples:
                    # Check if it has been previously saved
                    if treated in done:
                        continue
                    # Discard replicate and run data
                    treated_elements = treated.split('_')
                    run_value = treated_elements[-1]
                    del treated_elements[-1]
                    del treated_elements[2]
                    # Check if conditions match
                    if treated_elements == treated_conditions[j]:
                        # Save in group list
                        done.append(treated)
                        con = '_'.join(treated_elements)
                        groups_list[i].append(f'{run_value}/{con}')
        
    return groups_list


## MAIN PROGRAM

def main():
    '''
    Main program
    '''

    parser = argparse.ArgumentParser(prog='03-FilterByDepthAndRep.py', 
                                     description='''This program filters the
                                     libraries of a given species by selecting
                                     those that exceed the minimum thresholds
                                     of sequencing depth and number of 
                                     replicates specified by the user''',  
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-i', '--path-project', type=str, nargs=1,  
                        help='Absolute path to the directory in which the \
                              project libraries are located.')
    parser.add_argument('-o', '--path-out', type=str, nargs=1, 
                        help='Absolute path of the directory where the \
                            libraries (.fasta and .fastq) that pass the \
                            sequencing depth and number of replicates \
                            filters will be saved (within this directory \
                            the libraries will be stored in the subdirectory \
                            of the project they belong to).')
    parser.add_argument('-s', '--summary-dir',type=str, nargs=1,
                        help='Absolute path of the directory in which the \
                            summary files will be saved.')
    parser.add_argument('-m', '--metadata-dir', type=str, nargs=1,
                        help='Absolute path of the directory in which the \
                            metadata files are located.')
    parser.add_argument('-d', '--depth-threshold', type=int, nargs=1, 
                        help='Sequencing depth threshold. This argument must \
                            be a positive numerical value.')
    parser.add_argument('-r', '--rep-threshold', type=int, nargs=1,
                        help='Number of replicates threshold. This argument \
                            must be a positive numerical value.')
    parser.add_argument('-p', '--processes', type=int, nargs=1, default=1,
                        help='Number of processes (default: %(default)s)')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    
    ## 1. CHECK ARGUMENTS
    #######################################################################

    try:
        path_project = args.path_project[0]
        path_out = args.path_out[0]
        path_out_sum = args.summary_dir[0]
        path_metadata_dir = args.metadata_dir[0]
        depth_threshold = args.depth_threshold[0]
        rep_threshold = args.rep_threshold[0]
        num_process = args.processes[0]
    except:

        print('ERROR: You have inserted a wrong parameter or you are missing a parameter.')
        parser.print_help()
        sys.exit()
    
    # Create output summary directory
    os.system('mkdir -p ' + path_out_sum)
    
    # Get the species pathway and list their projects
    project = os.path.basename(path_project)
    species = os.path.basename(os.path.dirname(path_project))

    # Create output summary paths
    results_s_path = f'{path_out_sum}/{species}_{project}_sum_projects.tsv'
    summary_s_path = f'{path_out_sum}/{species}_{project}_sum_libraries.tsv'

    # Check if the project metadata file exists
    project_met_path = f'{path_metadata_dir}/{species}_m_{project}.tsv'
    if os.path.exists(project_met_path):

        ## 2. CHECK AND FILTER BY SEQUENCING DEPTH
        ###################################################################

        # Get project metadata table
        project_metadata = GetProjectMetadata(project_met_path)
        
        # Check the sequencing depth of the libraries
        depth_info_lib = GetLibrariesDepthPROCCESSING(path_project, project_metadata, depth_threshold, num_process)

        ## 3. CHECK AND FILTER BY NUMBER OF REPLICATES
        ###############################################################

        # Get a list of all the libraries in the project (run names)
        total_libraries = [run_name.replace('_tr.fastq.gz','') for run_name in depth_info_lib.keys()]

        # Get run names from file names
        lib_valid_names = [key.replace('_tr.fastq.gz','') for key, value in depth_info_lib.items() if value[1] == "valid"]

        # Obtain control-treatment groups from libraries selected for their sequencing depth.
        project_groups = GetProjectGroups(lib_valid_names, project_met_path)

        # Check if the control and treated conditions have the minimum number of replicates.
        filtered_libraries = []
        for group in project_groups:

            # Store the condition as a dictionary key and the runs
            # associated with that condition as a value (in a list).
            dic_conditions = {}
            for sample in group:
                sam_elements = sample.split("/") # SRR1848795/treated_drought_30%_T.0_leaves_27d
                run = sam_elements[0] # SRR1848795
                con = sam_elements[1] # treated_drought_30%_T.0_leaves_27d
                if con not in dic_conditions:
                    dic_conditions[con] = [run]
                else:
                    dic_conditions[con].append(run)
            
            # Iterate previously created dictionary
            valid_list_group = []
            valid_control = False
            valid_treated = False
            for condition in dic_conditions:

                # Count the number of libraries of the same condition
                num_samples = len(dic_conditions[condition])
                # If it is "control" and has more than rep_threshold replicates
                if 'control' in condition and num_samples >= rep_threshold:
                    valid_list_group += dic_conditions[condition]
                    valid_control = True
                # If it is "treated" and has more than rep_threshold replicates
                elif num_samples >= rep_threshold:
                    valid_list_group += dic_conditions[condition]
                    valid_treated = True

            # There must be at least "rep_threshold" control and treated replicates.
            if valid_control and valid_treated:
                filtered_libraries += valid_list_group

        # If valid libraries exist...
        if len(filtered_libraries) > 0:

            # Create project output directory
            path_out_project = f'{path_out}/fasta/{species}/{project}'
            path_out_project_fq = f'{path_out}/fastq/{species}/{project}'
            os.system(f'mkdir -p {path_out_project}')
            os.system(f'mkdir -p {path_out_project_fq}')

            # Save filtered libraries
            for lib_name in filtered_libraries:
                # Create fastq.gz and fasta paths
                path_fasta_lib = os.path.join(path_project, f'{lib_name}_tr.fasta')
                path_fastq_lib = os.path.join(path_project, f'{lib_name}_tr.fastq.gz')
                
                # Move fasta and fastq file to output directory
                os.system(f'mv {path_fasta_lib} {path_out_project}')
                os.system(f'cp {path_fastq_lib} {path_out_project_fq}')

                # Write the path of the selected files.
                with open(summary_s_path, 'a') as filtered:
                    # Get depth value from dictionary of selected files
                    fastq_gz_name = f'{lib_name}_tr.fastq.gz'
                    filtered.write(f'{path_fastq_lib}\t{depth_info_lib[fastq_gz_name][0]}\t{depth_info_lib[fastq_gz_name][1]}\tvalid\n')

        # Get a list of discarded libraries
        discarded_libraries = list(set(total_libraries) - set(filtered_libraries))

        # If discarded libraries exist...
        if len(discarded_libraries) > 0:
            # Save discarded libraries
            for lib_name in discarded_libraries:

                # Create fastq.gz and fasta paths
                path_fasta_lib = os.path.join(path_project, f'{lib_name}_tr.fasta')
                path_fastq_lib = os.path.join(path_project, f'{lib_name}_tr.fastq.gz')

                # Write the path of the discarded files.
                with open(summary_s_path, 'a') as excluded:
                    # Get depth value from dictionary of selected files
                    fastq_gz_name = f'{lib_name}_tr.fastq.gz'
                    excluded.write(f'{path_fastq_lib}\t{depth_info_lib[fastq_gz_name][0]}\t{depth_info_lib[fastq_gz_name][1]}\tnot-valid\n')
        
        # Save global results
        with open(results_s_path, 'a') as results:
            results.write(f'{species}\t{project}\t{str(len(filtered_libraries))}\t{str(len(discarded_libraries))}\n')

        # Delete discarded fasta files
        os.system(f'rm -f {path_project}/*.fasta')

    # If the project has no metadata it will be discarded
    else:
        # Save global results
        with open(results_s_path, 'a') as results:
            results.write(f'{species}\t{project}\tDiscarded\tDiscarded\n')
    
    print(f'{project} ({species}) done!')


## CALL THE MAIN PROGRAM

if __name__ == '__main__':
    '''
    Call to the main program
    '''
    main()
