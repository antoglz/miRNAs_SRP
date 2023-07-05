#!/bin/bash

#******************************************************************************
#  
#   07-Search_miRBase_by_taxonomy.sh
#
#   This program searches for sequences of a specific taxonomic group
#   in a .FASTA file with sequences downloaded from miRBase. Those
#   sequences that belong to the desired taxonomic group are saved in
#   a new fasta file whose name is specified by the user.
# 
#   Author: Antonio Gonzalez Sanchez
#   Date: 6/04/2023
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
    usage: 07-Search_miRBase_by_taxonomy.sh [options] ...
    options:
        -h                      help
        -d/--database           Path to the fasta file with the downloaded
                                sequences from the miRBase database.
        -o/--organism-table     Path to the table with information on
                                organisms from which miRNAs have been
                                identified (organisms.txt.gz file in miRBase).
        -t/--taxonomic-group    Taxonomic group whose sequences are desired
                                (e.g. Viridiplantae)
        -n/--output-name        Output file name (without extension)
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
   TEMP=$(getopt -o h::d:o:n:t: --long help::,database:,organism-table:,output-name:,taxonomic-group: -- "$@")

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
            -d|--database)
                database="$2"; shift 2 ;;
            -o|--organism-table)
                organism_table_path="$2"; shift 2 ;;
            -t|--taxonomic-group)
                tax="$2"; shift 2 ;;
            -n|--output-name)
                output_name="$2"; shift 2 ;;
            # -- meands the end of the arguments; drop this, and break out the while loop
            --) shift ; break ;;
            # If invalid options were passed...
            *) echo "Unexpected option: $1 - this should not happen."
                Usage ;;
      esac
   done
}


### MAIN
main() {

    # Get arguments
    ArgumentsManagement "$@"

    printf "Searching for $tax sequences...\n"

    # Read the input FASTA file
    while read line
    do
        # Find a header
        if [[ ${line:0:1} == '>' ]]
        then
            # Save the species name of the sequence in a variable
            seq_name=$(echo "$line" | awk '{print $3 " " $4}')

            # Search for the species name of the sequence in the organism table file.
            grep_result=$(grep -w "$seq_name" $organism_table_path)

            # If the species name of the sequence is found in the table...
            if [ $? -eq 0 ]
            then
                # If it contains the word 'Viridiplantae'
                if echo "$grep_result" | grep -q "$tax"
                then
                    # Write the header and sequence of the sequence to a new FASTA file.
                    echo "$line" >> $output_name.fa
                    read line
                    echo "$line" >> $output_name.fa
                fi
            fi
        fi
    done < $database

    printf "Done!\n"

}
main "$@"
