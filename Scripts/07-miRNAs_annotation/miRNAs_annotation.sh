#!/bin/bash

#******************************************************************************
#  
#   miRNAs_annotation.sh
#
#   This program annotates significant differentially expressed miRNAs
#   using the miRBase and PmiREN databases without considering a threshold
#   value of log2FC. To do this,it generates a fasta file with the sequences
#   and aligns them with the miRBase and PmiREN databases to identify the
#   mature miRNAs and also its precursors. Once annotated, the results of
#   both alignments are merged into a single five-column table: sequence,
#   Annotation with miRBase using only sequences of specific species,
#   Annotation with miRBase using the rest of the species, Annotation with
#   PmiREN using only sequences of specific species, Annotation with PmiREN
#   using the rest of the species. If any of the sequences has not obtained
#   results in the alignment with one of the two databases, it will be
#   represented as NULL. In addition, this program generates a summary file
#   with the number of sequences annotated for each of the two cases: miRNAs
#   and precursors.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 26/01/2023
#
#
#******************************************************************************



### FUNCTIONS

#################################################
#
#   This function contains the usage message
#
#################################################

usage() {
    
    help='''
    usage: miRNAs_annotation.sh [options] ...
    options:
        -h                      help
        -i/--input              Path to the directory where the species
                                directories are located. Tables with the
                                results of the differential expression
                                analysis can be found in each of these
                                directories.
        -o/--output             Path to the directory where the miRNA annotation
                                results will be stored.
        -s/--mirbase            Path to the directory where the downloaded
                                miRNAs and precursors (hairpin) files from
                                miRBase are located.
        -m/--pmiren             Path to the directory where the downloaded
                                miRNAs and precursors (hairpin) files from
                                Plant miRNA ENcyclopedia (PmiREN) are located.
        -d/--species-ids        Path to the file (table) that establishes the
                                correspondence between the three-letter code
                                and four-letter code used in this pipeline to
                                name the species.
        -r/--threads            Number of threads.
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

arguments_management() {

    ## 1. GET ARGUMENTS
    ############################################################################

    # Set a default value for optional arguments
    projects_list=''
    ea_table=''
    mww_pvalue=-1

    # Read the options
    TEMP=$(getopt -o h::i:o:m:e:s:t:l:a:p: --long help::,input:,output:,mirbase:,pmiren:,species-ids:,threads:,projects-list:,ea-table:,mww-pvalue: -- "$@")

    # Check if the arguments are valid
    VALID_ARGUMENTS=$?
    if [ "$VALID_ARGUMENTS" != 0 ];
    then
        usage
    fi

    # Convert options to script arguments
    eval set -- "$TEMP"

    # Extract options and their arguments into variables.
    while true
    do
        case "$1" in
            -h|--help)
                usage ;;
            -i|--input)
                path_in="$2"; shift 2 ;;
            -o|--output)
                path_out="$2"; shift 2 ;;
            -m|--mirbase)
                path_mirbase="$2"; shift 2 ;;
            -e|--pmiren)
                path_PmiREN="$2"; shift 2 ;;
            -s|--species-ids)
                path_ids_table="$2"; shift 2 ;;
            -t|--threads)
                threads="$2";
                # Check if -d is a positive number
                case "$threads" in
                    # Argument is not a positive number
                    ''|*[!0-9]*)
                        echo "Error: -t (--threads) must be a positive number"
                        usage ;;
                esac
                shift 2;;
            -l|--projects-list)
                # Check if a value for -l is provided
                if [ -n "$2" ] ; then
                    projects_valid_list="$2"
                fi
                shift 2 ;;
            -a|--ea-table)
                # Check if a value for -l is provided
                if [ -n "$2" ] ; then
                    ea_table="$2"
                fi
                shift 2 ;;
            -p|--mww-pvalue)
                # Check if a value for -p is provided
                if [ -n "$2" ]; then
                    mww_pvalue="$2"
                    # Check if -p is a positive number (including decimals)
                    if ! [[ "$mww_pvalue" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                        echo "Error: -p (--mww-pvalue) must be a positive number"
                        usage
                    fi
                fi
                shift 2 ;;
            # -- meands the end of the arguments; drop this, and break out the while loop
            --) shift ; break ;;
            # If invalid options were passed...
            *) echo "Unexpected option: $1 - this should not happen."
                usage ;;
        esac
    done


    ## 2. CHECK THE STRATEGY THAT WILL BE USED
    ############################################################################
    # Way 1 -> Use all the subprojects (Default)
    # Way 2 -> Use the subprojects from the list
    # Way 3 -> Use the subprojects that obtained a p-value (Exploratory Analysis)
    #          lower than 'mww_pvalue'

    # Default way
    way=1

    # Is mww_pvalue less than 0? (yes = 1, no = 0)
    result_pvalue=$(echo "$mww_pvalue < 0" | bc -l)

    # Check that the variable project_list is empty
    [[ -n $projects_valid_list ]] && way=2

    # Check that variable ea_table or variable mww_pvalue is not empty.
    [[ -n "$ea_table" || "$result_pvalue" -eq 0 ]] && way=3

    # Check that variables from way 2 and 3 have not been provided.
    [[ -n "$projects_valid_list" && ( "$result_pvalue" -eq 0 || -n "$ea_table" ) ]] &&
        { echo 'The optional arguments -l (--projects-list) and -p (--mww-pvalue) / -a (--ea-table) are mutually exclusive. Choose one option or none, but not both'; exit; }

    # Check if only one of the way 3 variables has been provided
    [[( -z "$ea_table" || "$result_pvalue" -eq 1 ) && "$way" -eq 3 ]] &&
        { echo "You cannot provide only one of the arguments ea-table and mww-pvalue. You must provide both or none"; exit; }
}


#####################################################
#
#   This function separates a database (.fasta)
#   into two different files, one with all the
#   sequences belonging to a specific species and
#   the other with the sequences of the rest of
#   the species in the database. To do this, it
#   receives the four-letter "id" of the species
#   in question, obtains the three-letter code of
#   that species establishing a correspondence
#   between both codes thanks to the file
#   'ids_table_file' and, finally, performs the
#   separation of the database file.
#
#   Arguments:
#       species id
#       Ids table file
#       Input file path
#       Output directory path
#
#####################################################

create_temporary_databases(){

    # Arguments
    local species_id="${1}"
    local ids_table_file="${2}"
    local path_file_in="${3}"
    local path_dir_out="${4}"

    # Other variables
    file_out_species=$path_dir_out/database_species.fasta
    file_out_others=$path_dir_out/database_rest.fasta
    path_dir_in=${path_file_in%/*}
    
    # Get miRBase and PmiREN ids
    id=$(grep $species_id $ids_table_file | awk -F"," '{print $2}')
    id_up=`echo ${id:0:1} | tr '[a-z]' '[A-Z]'`${id:1}
    
    # Replace uracils with thymines
    cat $path_file_in | cut -d " " -f 1 | awk 'BEGIN{RS=">"} {gsub("U", "T", $2); print ">"$1"\n"$2}' > $path_dir_in/temp_file.fasta
    tail -n +3 $path_dir_in/temp_file.fasta > $path_dir_in/temp2_file.fasta

    # Separate database by species and others.
    cat $path_dir_in/temp2_file.fasta | cut -d " " -f 1 | awk -v id="$id" \
                                              -v id_up="$id_up" \
                                              -v file_sp="$file_out_species" \
                                              -v file_o="$path_dir_in/temp_others_file.fasta" \
                                              'BEGIN{RS=">"} {if ($0 ~ id || $0 ~ id_up) {print ">"$1"\n"$2 > file_sp } else {print ">"$1"\n"$2 > file_o}}'

    # Remove the first two faulty lines from the 'others' file.
    tail -n +3 $path_dir_in/temp_others_file.fasta > $file_out_others

    # Remove temporary files
    rm $path_dir_in/*temp*
    
}


#####################################################
#
#   This function creates a fasta file from
#   the sequences located in the first column
#   of a results table obtained from a
#   differential expression analysis (DEA).
#   The name assigned to each sequence is
#   formed by the string 'sequence' followed
#   by a number that increments with each
#   sequence.
#
#   Arguments:
#       DEA results file path.
#       Output directory path.
#
#####################################################

table_to_fasta_and_tsv () {

    # Arguments
    local path_file="${1}"
    local path_out="${2}"

    # Get filename
    file=${path_file##*/}
    filename=${file%.csv}
    seqnum=1

    # Read file
    while IFS=, read seq other
    do
        # Discard header
        if [[ $seq != "seq" ]]
        then
            # Create fasta file with sequences
            echo ">sequence"$seqnum >> $path_out/$filename.fasta
            echo $seq >> $path_out/$filename.fasta

            # Create tsv file with sequences
            printf "sequence"$seqnum"\t"$seq"\n" >> $path_out/$filename.tsv
            let "seqnum++" 
        fi
    done < $path_file

}


#####################################################
#
#   This function joins two TSV tables by
#   INNER JOIN or FULL OUTER JOIN, generating
#   an output file with the same format. It
#   receives the absolute paths of both TSV
#   files, the key column of each table to
#   perform the join, the columns of each
#   table to be saved in the output file,
#   the absolute path of the output file and
#   an argument specifying whether the
#   procedure will be done by INNER JOIN or
#   FULL OUTER JOIN.
#
#   Arguments:
#       Absolute path of TSV file 1
#       Absolute path of TSV file 1
#       Key column of TSV file 1
#       Key column of TSV file 2
#       Output file path
#       String specifying whether the join
#           will be made by INNER JOIN
#           (FALSE) or FULL OUTER JOIN (TRUE)
#           (Optional).
#
#####################################################

merge_tsv_files () {

     # Arguments
    local path_file_1="${1}"
    local path_file_2="${2}"
    local key_file_1="${3}"
    local key_file_2="${4}"
    local output_cols="${5}"
    local output_file="${6}"
    local outer_join="${7:-FALSE}"

    # Files name
    filename1_path=${path_file_1%.tsv}
    filename2_path=${path_file_2%.tsv}

    # Sort files to be joined
    LANG=en_EN sort -k $key_file_1 -t$'\t' $path_file_1 -o $filename1_path"_sort.tsv"
    LANG=en_EN sort -k $key_file_2 -t$'\t' $path_file_2 -o $filename2_path"_sort.tsv"
    
    # Inner join
    if [ $outer_join == "FALSE" ]
    then
        LANG=en_EN join -1 $key_file_1 -2 $key_file_2 -t$'\t' \
            -o $output_cols -e "NULL" \
            $filename1_path"_sort.tsv" \
            $filename2_path"_sort.tsv"  > $output_file
    # Full outer join
    elif [ $outer_join == "TRUE" ]
    then
        LANG=en_EN join -1 $key_file_1 -2 $key_file_2 -t$'\t' \
            -o $output_cols -e "NULL" \
            -a 1 -a 2 \
            $filename1_path"_sort.tsv" \
            $filename2_path"_sort.tsv"  > $output_file
    else
        echo "Error: Parameter incorrectly entered in merge_tsv_files function!"
        echo "Please specify whether the union of TSV files will be done by INNER JOIN (FALSE) or FULL OUTER JOIN (TRUE)."
    fi

}


#####################################################
#
#   This function creates an annotation table
#   from 4 SAM files obtained from the
#   alignment of a fasta file with the miRBase
#   and PmiREN databases (sequences belonging
#   to the species and sequences belonging to
#   other species separately). This table
#   consists of 5 columns: seq, miRBase_
#   species, miRBase_rest, PmiREN_species,
#   PmiREN_rest. The first column contains
#   the nucleotide sequence and the last 4
#   columns the name of each sequence according
#   to the miRBase (species and other species)
#   and PmiREN (species and other species)
#   databases, respectively.
#
#   Arguments:
#       SAM file obtained from miRBase
#           alignment (species)
#       SAM file obtained from miRBase
#           alignment (other species)
#       SAM file obtained from PmiREN
#           alignment (species)
#       SAM file obtained from PmiREN
#           alignment (other species)
#       Output file name
#       Output directory path.
#
#####################################################

get_miRNAs_annotation () {

    # Arguments
    local sam_file_mirbase_species="${1}"
    local sam_file_mirbase_rest="${2}"
    local sam_file_PmiREN_species="${3}"
    local sam_file_PmiREN_rest="${4}"
    local path_tsv_seq="${5}"
    local out_name="${6}"
    local path_out="${7}"

    # Get annotation tables
    [[ -f $sam_file_mirbase_species ]] &&
        cat $sam_file_mirbase_species | cut -f1,3 | awk -F "\t" '{if (substr($1,1,1) != "@"){print}}' > $path_out/$out_name"_mirbase_species_temp.tsv" ||
        touch $path_out/$out_name"_mirbase_species_temp.tsv"
    [[ -f $sam_file_mirbase_rest ]] && 
        cat $sam_file_mirbase_rest | cut -f1,3 | awk -F "\t" '{if (substr($1,1,1) != "@"){print}}' > $path_out/$out_name"_mirbase_rest_temp.tsv" ||
        touch $path_out/$out_name"_mirbase_rest_temp.tsv" 
    [[ -f $sam_file_PmiREN_species ]] &&
        cat $sam_file_PmiREN_species | cut -f1,3 | awk -F "\t" '{if (substr($1,1,1) != "@"){print}}' > $path_out/$out_name"_PmiREN_species_temp.tsv" ||
        touch $path_out/$out_name"_PmiREN_species_temp.tsv"
    [[ -f $sam_file_PmiREN_rest ]] &&
        cat $sam_file_PmiREN_rest | cut -f1,3 | awk -F "\t" '{if (substr($1,1,1) != "@"){print}}' > $path_out/$out_name"_PmiREN_rest_temp.tsv" ||
        touch $path_out/$out_name"_PmiREN_rest_temp.tsv"

    # Join nucleotide sequence with identifier
    merge_tsv_files $path_out/$out_name"_mirbase_species_temp.tsv" $path_tsv_seq 1 1 1.2,2.2 $path_out/$out_name"_mirbase_species_annot.tsv"
    merge_tsv_files $path_out/$out_name"_mirbase_rest_temp.tsv" $path_tsv_seq 1 1 1.2,2.2 $path_out/$out_name"_mirbase_rest_annot.tsv"
    merge_tsv_files $path_out/$out_name"_PmiREN_species_temp.tsv" $path_tsv_seq 1 1 1.2,2.2 $path_out/$out_name"_PmiREN_species_annot.tsv"
    merge_tsv_files $path_out/$out_name"_PmiREN_rest_temp.tsv" $path_tsv_seq 1 1 1.2,2.2 $path_out/$out_name"_PmiREN_rest_annot.tsv"
    
    # Create tsv annotation table
    merge_tsv_files $path_out/$out_name"_mirbase_species_annot.tsv" $path_out/$out_name"_mirbase_rest_annot.tsv" 2 2 0,1.1,2.1 $path_out/$out_name"_temp_1.tsv" TRUE
    merge_tsv_files $path_out/$out_name"_temp_1.tsv" $path_out/$out_name"_PmiREN_species_annot.tsv" 1 2 0,1.2,1.3,2.1 $path_out/$out_name"_temp_2.tsv" TRUE
    merge_tsv_files $path_out/$out_name"_temp_2.tsv" $path_out/$out_name"_PmiREN_rest_annot.tsv" 1 2 0,1.2,1.3,1.4,2.1 $path_out/$out_name"_temp_3.tsv" TRUE

    # Convert tsv to csv
    sort -k 1b,1 -t$'\t' $path_out/$out_name"_temp_3.tsv" -o $path_out/$out_name"_temp_3_sort.tsv"
    echo "seq,miRBase_species,miRBase_others,PmiREN_species,PmiREN_others" > $path_out/$out_name"_annot.csv"
    cat $path_out/$out_name"_temp_3_sort.tsv" | sed $'s/\t/,/g' >> $path_out/$out_name"_annot.csv"
  
    # Remove temporal files
    rm $path_out/*temp*
    rm $path_out/*PmiREN*
    rm $path_out/*mirbase*
}


#####################################################
#
#   This function filters annotated sequences,
#   selecting only those that have been annotated
#   at least once in both the miRBase and PmiREN
#   databases. It processes the CSV file
#   corresponding to the annotation table
#   containing the following fields: seq,
#   miRBase_species, miRBase_others, PmiREN_
#   species, and PmiREN_others. It checks each
#   line in the input file to determine if the
#   miRBase fields (miRBase_species and miRBase_
#   others) or the PmiREN fields (PmiREN_species
#   and PmiREN_others) are not equal to "NULL".
#   If at least one field in both miRBase and
#   PmiREN is not "NULL", the functions writes
#   the line to a new file (path_file_out). This
#   function effectively filters and saves lines
#   from the input CSV file based on specified
#   criteria.
#
#   Arguments:
#       Input Annotation File Path.
#       Output Path for Filtered Annotation File.
#
#####################################################

filter_annotated_sequences () {

    # Arguments
    local path_file_in="${1}"
    local path_file_out="${2}"

    # Read the original table and iterate through each line
    while IFS=, read -r seq miRBase_species miRBase_others PmiREN_species PmiREN_others
    do  
        # Check if the sequence is annotated at least once in miRBase and PmiREN
        if [[ ( "$miRBase_species" != "NULL" || "$miRBase_others" != "NULL" ) \
            && ( "$PmiREN_species" != "NULL" || "$PmiREN_others" != "NULL" ) ]]; then
            
            # Print the line that meets the criteria and save it to the output file
            echo "$seq,$miRBase_species,$miRBase_others,$PmiREN_species,$PmiREN_others" >> "$path_file_out"
        fi
    done < "$path_file_in"

}


#####################################################
#
#   This function takes the absolute path of the
#   summary table obtained from the Diff_exp_
#   analysis.r program and filters projects based
#   on a p-value threshold (MWW). It returns a
#   list of valid projects.
#
#   Arguments:
#       Input ea_summary.csv file path.
#       MWW p-value threshold
#
#####################################################

filter_subprojects_by_mww_pvalue (){
    
    # Arguments
    local path_file_in="${1}"
    local mww_pvalue="${2}"

    # Filter sequences by pvalue (MWW)
    filtered_projects=$(tail -n +2 $path_file_in | awk -F ',' -v var=$mww_pvalue '{if ($NF <= var) print $2}')
    
    # Return
    echo $filtered_projects

}


### MAIN
main () {

    # Get arguments
    arguments_management "$@"

    # Way 3: Get the valid subprojects from the provided list
    [ "$way" -eq 3 ] && valid_subprojects=$(filter_subprojects_by_mww_pvalue "$ea_table" "$mww_pvalue")

    # Get directory name
    data=$(basename "$path_in")

    # Get data suffix (raw or sig)
    suffix=$(echo $data | awk -F '_' '{print $2}')

    # Iterate reference list
    references_list=$(ls $path_mirbase)
    for reference in $references_list
    do
        # Create output directory and summary file path
        reference_name=${reference%.fa} # hairpin.fa or mature.fa
        mkdir -p $path_out/miRNA_annotation_$suffix/$reference_name
        path_out_summary=$path_out/miRNA_annotation_$suffix/$reference_name/summary.csv

        # Iterate species
        species_list=$(ls $path_in)
        for species in $species_list
        do
            # Create projects_path and list projects
            projects_path=$path_in/$species
            projects_list=$(ls $projects_path)

            # Create temporary databases (species and others)
            create_temporary_databases $species $path_ids_table $path_mirbase/$reference $path_mirbase
            create_temporary_databases $species $path_ids_table $path_PmiREN/$reference $path_PmiREN

            # Create Index directories
            mkdir -p $path_mirbase/mirbase_species_idx
            mkdir -p $path_mirbase/mirbase_rest_idx
            mkdir -p $path_PmiREN/PmiREN_species_idx
            mkdir -p $path_PmiREN/PmiREN_rest_idx

            # Index databases
            bowtie2-build --threads $threads $path_mirbase/database_species.fasta $path_mirbase/mirbase_species_idx/mirbase_species > /dev/null 2>&1
            bowtie2-build --threads $threads $path_mirbase/database_rest.fasta $path_mirbase/mirbase_rest_idx/mirbase_rest > /dev/null 2>&1
            bowtie2-build --threads $threads $path_PmiREN/database_species.fasta $path_PmiREN/PmiREN_species_idx/PmiREN_species > /dev/null 2>&1
            bowtie2-build --threads $threads $path_PmiREN/database_rest.fasta $path_PmiREN/PmiREN_rest_idx/PmiREN_rest > /dev/null 2>&1
            
            # Iterate projects
            for project in $projects_list
            do

                # Create projects_path and list projects
                files_path=$projects_path/$project
                files_list=$(ls $files_path)

                # Iterate files
                for file in $files_list
                do

                    # Output name
                    out_name=${file%_dea*}

                    # Way 2: Filtering using subprojects list
                    if [ "$way" -eq 2 ]
                    then
                        # Check if subproject is valid
                        check_subproject_valid=$( grep $file $projects_valid_list )

                        # If it is not valid, pass to the next subproject
                        [[ -z $check_subproject_valid ]] && valid="false" || valid="true"

                        # Save results in summary file and pass to the next subproject.
                        if [[ $valid == "false" ]]
                        then
                            echo $species","$out_name",NA,NA" >> $path_out_summary
                            continue
                        fi

                    # Way 3: Filtering using MWW p-value
                    elif [ "$way" -eq 3 ]
                    then
                        # Delete suffix
                        file_suffix_removed=$(echo "$file" | sed 's/_dea_sig\.csv$//')

                        # If it is not valid, pass to the next subproject
                        [[ $valid_subprojects != *$file_suffix_removed* ]] && valid="false" || valid="true"
                            
                        # Save results in summary file and pass to the next subproject.
                        if [[ $valid == "true" ]]
                        then
                            echo $species","$out_name",NA,NA" >> $path_out_summary
                            continue
                        fi
                    fi

                    # Output paths
                    path_fasta_files=$path_out/miRNA_annotation_$suffix/$reference_name/01-Significant_seq_fasta/$species/$project
                    path_sam_files=$path_out/miRNA_annotation_$suffix/$reference_name/02-Alignment_miRBase-PmiREN_results/$species/$project
                    path_annotation_files=$path_out/miRNA_annotation_$suffix/$reference_name/03-Annotation_tables/$species/$project
                    path_annotation_files_filt=$path_out/miRNA_annotation_$suffix/$reference_name/04-Annotation_tables_filtered/$species/$project

                    # Create output directories
                    mkdir -p $path_fasta_files
                    mkdir -p $path_sam_files
                    mkdir -p $path_annotation_files
                    mkdir -p $path_annotation_files_filt

                    printf "\n########################### File: $file ($species) ########################### \n"
        

                    # Remove header and create temporary file
                    tail -n +2 $files_path/$file > $files_path/temp_file.csv
                    
                    # Check if the file does not contain differentially expressed sRNAs (file empty)
                    if [ -s $files_path/temp_file.csv ]
                    then

                        # Create Fasta
                        table_to_fasta_and_tsv $files_path/$file $path_fasta_files

                        # Check if fasta file exist
                        if test -f "$path_fasta_files/$out_name""_dea_"$suffix".fasta"
                        then

                            # Alignment with miRBase (species database)
                            printf "Bowtie2 alignment with miRBase (species database)...\n"
                            results_mirbase_sp=$(bowtie2 -f --no-unal -p $threads \
                                                         -x $path_mirbase/mirbase_species_idx/mirbase_species \
                                                         -U $path_fasta_files/$out_name"_dea_"$suffix".fasta" \
                                                         -S $path_sam_files/$out_name"_mirbase_species.sam" 2>&1)
                            
                            # Check if an alignment error occurred due to the non-existence of the database.
                            [ $? -ne 0 ] && echo "Species $species not found in the miRBase database!" || printf "Done!\n"

                            # Alignment with miRBase (Rest database)
                            printf "Bowtie2 alignment with miRBase (others database)...\n"
                            results_mirbase_others=$(bowtie2 -f --no-unal -p $threads \
                                                             -x $path_mirbase/mirbase_rest_idx/mirbase_rest \
                                                             -U $path_fasta_files/$out_name"_dea_"$suffix".fasta" \
                                                             -S $path_sam_files/$out_name"_mirbase_rest.sam" 2>&1)
                                                    
                            # Check if an alignment error occurred due to the non-existence of the database.
                            [ $? -ne 0 ] && echo "Species $species not found in the miRBase database!" || printf "Done!\n"

                            # Alignment with PmiREN (species database)
                            printf "Bowtie2 alignment with PmiREN (species database)...\n"
                            results_pmiren_sp=$(bowtie2 -f --no-unal -p $threads \
                                                        -x $path_PmiREN/PmiREN_species_idx/PmiREN_species \
                                                        -U $path_fasta_files/$out_name"_dea_"$suffix".fasta" \
                                                        -S $path_sam_files/$out_name"_PmiREN_species.sam" 2>&1)

                            # Check if an alignment error occurred due to the non-existence of the database.
                            [ $? -ne 0 ] && echo "Species $species not found in the PmiREN database!" || printf "Done!\n"
                            
                            # Alignment with PmiREN (Rest database)
                            printf "Bowtie2 alignment with PmiREN (others database)...\n"
                           results_pmiren_others=$(bowtie2 -f --no-unal -p $threads \
                                                           -x $path_PmiREN/PmiREN_rest_idx/PmiREN_rest \
                                                           -U $path_fasta_files/$out_name"_dea_"$suffix".fasta" \
                                                           -S $path_sam_files/$out_name"_PmiREN_rest.sam" 2>&1)
                            
                            # Check if an alignment error occurred due to the non-existence of the database.
                            [ $? -ne 0 ] && echo "Species $species not found in the PmiREN database!" || printf "Done!\n"
                            
                            # Create annotation table from sam files
                            printf "Creating annotation table...\n"
                            get_miRNAs_annotation $path_sam_files/$out_name"_mirbase_species.sam" \
                                                $path_sam_files/$out_name"_mirbase_rest.sam" \
                                                $path_sam_files/$out_name"_PmiREN_species.sam" \
                                                $path_sam_files/$out_name"_PmiREN_rest.sam" \
                                                $path_fasta_files/$out_name"_dea_"$suffix".tsv" \
                                                $out_name $path_annotation_files
                            printf "Done!\n"
                            
                            printf "Creating filtered annotation tables....\n"
                            filter_annotated_sequences $path_annotation_files/$out_name"_annot.csv" \
                                                     $path_annotation_files_filt/$out_name"_annot_filt.csv"
                            printf "Done!\n"


                            # Count the number of annotated sequences and how many of them have been selected.
                            num_miRNAs=$(tail -n +2 $path_annotation_files/$out_name"_annot.csv" | wc -l)
                            num_miRNAs_filtered=$(tail -n +2 $path_annotation_files_filt/$out_name"_annot_filt.csv" | wc -l)

                            # Save it in summary file
                            echo $species","$out_name","$num_miRNAs","$num_miRNAs_filtered >> $path_out_summary
                            
                            
                        else
                            printf "NOTE: $out_name""dea_"$suffix".fasta does not exist\n"
                            echo $species","$out_name",0,0" >> $path_out_summary
                        fi

                    # The file has no differentially expressed sRNAs. 
                    else
                        printf "NOTE: the file has no differentially expressed sRNAs (File is empty)\n"
                        echo $species","$out_name",0,0" >> $path_out_summary
                    fi
                    # Delete temporary file
                    rm $files_path/temp_file.csv
                done
                # end file loop
            done
            # end project loop

            # Remove database files
            rm -r $path_mirbase/*species*
            rm -r $path_mirbase/*rest*
            rm -r $path_PmiREN/*species*
            rm -r $path_PmiREN/*rest*
        done
        # end species loop
    done
    # end reference loop

    # Create final summary file
    path_miRNAs_summary=$path_out/miRNA_annotation_$suffix/mature/summary.csv
    path_precursor_summary=$path_out/miRNA_annotation_$suffix/hairpin/summary.csv

    # Sort files to be joined
    LANG=en_EN sort -k 2 -t ',' $path_miRNAs_summary -o $path_out/miRNA_annotation_$suffix/mature/summary"_sort.csv"
    LANG=en_EN sort -k 2 -t ',' $path_precursor_summary -o $path_out/miRNA_annotation_$suffix/hairpin/summary"_sort.csv" 

    # Join summary files
    LANG=en_EN join -1 2 -2 2 -t ',' \
            -o 1.1,1.2,1.3,1.4,2.3,2.4  \
            $path_out/miRNA_annotation_$suffix/mature/summary"_sort.csv" \
            $path_out/miRNA_annotation_$suffix/hairpin/summary"_sort.csv"  > $path_out/miRNA_annotation_$suffix/summary.csv
    sed -i '1 i\Species,Experiment,Annot_miRNAs,Annot_miRNAs_filtered,Annot_precursor,Annot_precursor_filtered' $path_out/miRNA_annotation_$suffix/summary.csv

    # Delete temporary files
    rm $path_out/miRNA_annotation_$suffix/mature/*summary*
    rm $path_out/miRNA_annotation_$suffix/hairpin/*summary*

}
main "$@"
