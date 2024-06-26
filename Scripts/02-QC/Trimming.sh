#! /bin/bash

#******************************************************************************
#  
#   02-Trimming.sh
#
#   This program is used to trim the sequences using Illumina adapters.
#   For this purpose, the Fastp program is utilized.
#
#   Author: Antonio Gonzalez Sanchez and Pascual Villalba
#   Date: 19/10/2023
#   Version: 2 
#
#******************************************************************************


### FUNCTIONS

#####################################################
#
#   This function manages the input of arguments
#   to the program.
#
#####################################################

ArgumentsManagement() {
   
   # Read the options
   TEMP=$(getopt -o h::p:o:a:n:x: --long help::,project:,output:,adapters:,minlength:,maxlength: -- "$@")

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
                path_in="$2"; shift 2 ;;
            -o|--output)
                path_out="$2"; shift 2 ;;
            -a|--adapters)
                path_adapters="$2"; shift 2 ;;
            -n|--minlength)
                min_len="$2"; shift 2 ;;
            -x|--maxlength)
                max_len="$2"; shift 2 ;;
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
#   This function trims the sequences from a
#   specific library using the fastp program.
#
#   Arguments:
#       Fastq.gz file path
#       Output directory path
#       Path to the file with Illumina adapters.
#
#####################################################

FastpTrimming(){

    # Arguments
    local path_in="${1}"
    local path_out="${2}"
    local path_adapters="${3}"
    local min_len="${4}"
    local max_len="${5}"

    # Get the input directory path
    path_in_dir=$(dirname "$path_in")

    # Get the file name and SRR
    file=$(basename "$path_in")
    srr=$(echo "$file" | awk -F '[_.]' '{print $1}')

    # Trimming
    ## Paired end
    if [[ "$file" == *_1.fastq.gz ]];
    then
        fastp --adapter_fasta $path_adapters \
            -i $path_in_dir/$srr"_1.fastq.gz" -I $path_in_dir/$srr"_2.fastq.gz" \
            -o $path_out/$srr"_tr_1P.fastq.gz" -O $path_out/$srr"_tr_2P.fastq.gz" \
            --correction --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
            --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 \
            --length_required $min_len --trim_poly_x --poly_x_min_len 10 \
            --length_limit $max_len --n_base_limit 0

    ## Single end
    elif [[ "$file" != *_2.fastq.gz ]];
    then
        fastp --adapter_fasta $path_adapters \
            -i $path_in \
            -o $path_out/$srr"_tr.fastq.gz" \
            --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
            --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 \
            --length_required $min_len --trim_poly_x --poly_x_min_len 10 \
            --length_limit $max_len --n_base_limit 0
    fi
}


### MAIN
main () {

    # Get arguments
    ArgumentsManagement "$@"

    # Get species and project names
    project=$(basename "$path_in")
    species=$(basename "$(dirname "$path_in")")

    # Create output directories
    path_out_project=$path_out/$species/$project
    mkdir -p $path_out_project

    # List and iterate fastq files of the project
    fastq_list=$( ls -d1 $path_in/*fastq.gz )
    for fastq in $fastq_list
    do
        # Trimming
        FastpTrimming $fastq $path_out_project $path_adapters $min_len $max_len
    done
}
main "$@"
