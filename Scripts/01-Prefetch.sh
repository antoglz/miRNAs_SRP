#!/bin/bash

#******************************************************************************
#  
#   01-Prefetch.sh
#
#   This program downloads the libraries associated with a specific SRA
#   project using prefetch and fasterq-dump. Additionally, it performs
#   quality control on the library using FastQC.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 05/04/2022
#   Version: 1.0 
#
#******************************************************************************


### FUNCTIONS

#################################################
#
#   This function contains the usage message
#
#################################################

Usage() {
    
    help='''
    usage: 01-Prefetch.sh [options] ...
    options:
        -h                      help
        -i/--input              Path to the ".txt" file containing the SRA (Run)
                                identifiers of a specific project.
        -o/--output             Path to the directory in which the libraries
                                and FastQC results will be stored.
    '''
    echo "$help"
    exit 0
}


#####################################################
#
#   This function manages the input of arguments
#   to the program.
#
#####################################################

ArgumentsManagement() {
   
   # Read the options
   TEMP=$(getopt -o h::p:o: --long help::,project:,output: -- "$@")

   # Check if the arguments are valid
   VALID_ARGUMENTS=$?
   if [ "$VALID_ARGUMENTS" != 0 ];
   then
      Usage
   fi

   # Convert options to script arguments
   eval set -- "$TEMP"

   # Extract options and their arguments into variables.
   while true
   do
        case "$1" in
            -h|--help)
                Usage ;;
            -p|--project)
                project_path="$2"; shift 2 ;;
            -o|--output)
                path_out="$2"; shift 2 ;;
            # -- meands the end of the arguments; drop this, and break out the while loop
            --) shift ; break ;;
            # If invalid options were passed...
            *) echo "Unexpected option: $1 - this should not happen."
                Usage ;;
      esac
   done
}


#####################################################
#
#   This program downloads the library associated
#   with a specific SRA run using prefetch and
#   fasterq-dump. Additionally, it performs
#   quality control on the library using FastQC.
#   The function takes the run identifier and
#   the output directory path as arguments, where
#   the libraries and FastQC results will be stored.
#
#   Arguments:
#       SRR path
#       Output path
#
#####################################################

DownloadAndFastqc () {

    # Arguments
    local srr="${1}"
    local path_out="${2}"

    # Download file
    prefetch $srr --verify yes --check-all
    mv ./$srr/$srr.sra .
    rm -r $srr

    # Extract the fastq file (single) or files (paired)
    fasterq-dump --split-files -O $path_out $srr.sra
    rm $srr.sra

    # Obtain fastqc data
    if [ -e $path_out/$srr"_1.fastq" ] && [ -e $path_out/$srr"_2.fastq" ];
    then
        # Compress to .gz.
        gzip $path_out/$srr"_1.fastq"
        gzip $path_out/$srr"_2.fastq"

        # Fastqc.
        fastqc $path_out/$srr"_1.fastq.gz"
        fastqc $path_out/$srr"_2.fastq.gz"
        mv $path_out/$srr"_1_fastqc.zip" $path_out/fastqc_data
        mv $path_out/$srr"_1_fastqc.html" $path_out/fastqc_data
        mv $path_out/$srr"_2_fastqc.zip" $path_out/fastqc_data
        mv $path_out/$srr"_2_fastqc.html" $path_out/fastqc_data
    else
        # Compress to .gz.
        gzip $path_out/$srr".fastq"

        # Fastqc.
        fastqc $path_out/$srr".fastq.gz"
        mv $path_out/$srr"_fastqc.zip" $path_out/fastqc_data
        mv $path_out/$srr"_fastqc.html" $path_out/fastqc_data
    fi
}


### MAIN
main () {
    
    # Get arguments
    ArgumentsManagement "$@"

    # Get file name and path to input directory
    project_file_name=$(basename "$project_path")
    path_in=$(dirname "$project_path")

    ## Name of directory and subdirectories
    species=$(echo "$project_file_name" | cut -d "_" -f 1)
    project=$(echo "$project_file_name" | cut -d "_" -f 2 | cut -d "." -f 1)

    ## Create directory and subdirectories
    path_out_project=$path_out/$species/$project
    mkdir -p $path_out_project/fastqc_data
    
    ## Download SRR and obtain fastqc data
    list_of_fields="$(cat $path_in/$project_file_name)"
    for SRR in $list_of_fields;
    do
        DownloadAndFastqc $SRR $path_out_project
    done
}
main "$@"