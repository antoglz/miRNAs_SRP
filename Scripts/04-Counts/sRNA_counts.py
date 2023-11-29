#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#******************************************************************************
#  
#   04-Counts_and_sRNADatabase.py
#
#   This program generates the absolute counts and Reads Per Million (RPM)
#   tables for a specific project using the trimmed and filtered libraries,
#   also calculating the averages of both types of counts in the different
#   replicates of each condition for the sequences analysed. To do this,
#   a series of procedures are carried out:
#   
#   1. Filter the libraries by RNAcentral (Optional)
#
#   Align the sequences of each library with the rRNA, tRNA, snRNA and
#   snoRNA sequences contained in the RNAcentral database for
#   identification and elimination.
#
#   2. Create the absolute counts and RPM tables for each library.
#
#   Read filtered libraries and create absolute counts and RPM tables
#   in CSV format. These tables have two columns: seq and counts (absolute
#   counts)/RPM (Reads per Million).
#   
#               RPM = absolute count * 1000000 / size library.
#  
#   3. Join all the RPM tables and all the absolute counts tables in
#      two tables.
#
#   Once the absolute counts and RPM have been calculated for the sequences
#   of each library, a table is created for each type of count by joining
#   the results of the libraries belonging to the same project.
#
#
#   Authors: Antonio Gonzalez Sanchez, Pascual Villalba Bermell (pvbermell)
#   Date: 20/09/2023
#   Version: 2.0
#
#******************************************************************************


## IMPORT MODULES

import argparse
import csv
import itertools
import multiprocessing
import numpy as np
from numpy import genfromtxt
import os
import pandas as pd
import sqlite3
from sqlite3 import Error
import sys
from time import time
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


## FUNCTIONS

### 1. SQLITE CONNECT, INSERT AND JOIN FUNCTIONS
def connect_to_database (database_name: str):
    '''
    This function establishes a connection to an SQlite database, creating it
    if it does not exist.

    Parameters
    ----------
    database_name : str
        Absolute database path

    Returns
    -------
    Connection
        Conexion
    Cursor
        Cursor
    '''
    try:
        ## Connect to SQlite
        sqliteConnection = sqlite3.connect(database_name)
        cursor = sqliteConnection.cursor()
        
    except Error as error:
        print('ERROR. Unable to connect to the database')
        print(f'ERROR:\n{error}\n')
        sys.exit()

    else:
        return sqliteConnection, cursor

def insert_to_database (database_name: str, table_name: str, data_path: str,
                        type_data: str='sequences') -> None:
    '''
    This function inserts data into a specified SQLite database table.

    Parameters
    ----------
    database_name : str
        SQLite database name (filename)
    table_name : str
        Name of the table in which the data will be inserted
    data_path : str
        Absolute path of the file that contains the data to insert
    type_data : str
        Type of data to be inserted. In this program you can insert a library
        of sequences in .tsv format (type_data="sequences"), a table of
        absolute counts (type_data="counts_abs") or a table of rpm
        ("type_data="counts_rpm"). By default "sequences.
    '''
    # Connect to database
    sqliteConnection, cursor = connect_to_database(database_name)

    # Select separator and type of data
    if type_data == 'sequences':
        sep ='\t'
        query_section = '(seq TEXT, header TEXT);'
    elif type_data == 'counts_abs':
        sep =','
        query_section = '(seq TEXT, counts INT);'
    elif type_data == 'counts_rpm':
        sep =','
        query_section = '(seq TEXT, RPM REAL);'
    
    # Read data to insert
    chunksize=1000000
    data_to_insert = pd.read_csv(data_path, sep=sep, chunksize=chunksize, low_memory=False)

    try:
        # Pragma adjust
        cursor.execute('PRAGMA synchronous = OFF')
        cursor.execute('PRAGMA journal_mode = OFF')

        # Create table
        cursor.execute(f'CREATE TABLE IF NOT EXISTS {table_name}{query_section}')

        # Insert chunks in table
        for chunk in data_to_insert:
            # Insert chunk data into table
            chunk.to_sql(name=table_name, con=sqliteConnection, if_exists='append', index=False)       
        print(f'Data inserted in {table_name} correctly!')
        sys.stdout.flush()

        # Create table index to increase query speed
        print(f'Creating index of {table_name}')
        sys.stdout.flush()
        cursor.execute(f'CREATE INDEX idx_{table_name} ON {table_name} (seq);')
        print('DONE!\n')
        sys.stdout.flush()
        exit = False

    except Error as error:
        print('ERROR: Something went wrong during the data insertion.')
        print(f'ERROR:\n{error}\n')
        exit = True
    
    finally:
        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()
        
        # If it fails, exit the program
        if exit:
            sys.exit()

def merge_counts_tables_processing (list_of_tuples: list) -> None:
    '''
    This function is used to parallelize the merge_counts_tables function,
    which generates and executes the necessary SQLite queries to join in
    the same table the different replicates of the same condition or multiple
    tables of different conditions that contain the data of all their
    replicates.

    Parameters
    ----------
    list_of_tuples : list
        List of two tuples. The first tuple contains the information associated
        with the absolute counts and the second one with the RPM. Each tuple
        contains 4 elements, which serve as arguments to the merge_counts_tables()
        function. The first of these arguments is a list with the names of the
        replicates of each condition (t_1_r_1, t_1_r_2, t_2_r_1, t_2_r_2), the
        second one is a string that specifies the type of data we are working
        with (it can be "counts" or "RPM"), the third one is a string that
        specifies the type of tables we want to join (it can be "replicates" or
        "conditions") and the last one is a string that specifies the method by
        which replicates of the same condition will join if the third argument
        is "replicates". The last element is optional.
    '''

    # Two processes will be used, one for absolute counts and one for RPM.
    num_process = 2

    # Processes list
    processes = []

    for i in range(num_process):
        # Get arguments stored in the tuple
        database_name = list_of_tuples[i][0]
        list_tables = list_of_tuples[i][1]
        type_data = list_of_tuples[i][2]
        type_tables = list_of_tuples[i][3]
        mode = list_of_tuples[i][4]
        
        # Write process
        processes.append(multiprocessing.Process(
                        target = merge_counts_tables,
                        args=(database_name, list_tables, type_data,
                              type_tables, mode)))
        
        # Execute process
        processes[i].start()
        print(f'Process {i + 1} launched.')
    
    # Wait for all processes to finish
    for p in processes:
        p.join()

def merge_counts_tables (database_name: str, data_in: list, type_data: str,
                         type_tables: str, mode: str='inner') -> None:
    '''
    This function generates and executes the necessary SQLite queries to join in
    the same table the different replicates of the same condition or multiple
    tables of different conditions that contain the data of all their replicates.
    The first of the two processes can be performed using the INNER JOIN (PCA)
    or FULL OUTER JOIN (Differential Expression Analysis) methods. However, the
    second one is done by default using the FULL OUTER JOIN method.

    Parameters
    ----------
    database_name : str
        Absolute path to the database where the tables to be joined are
        located.
    data_in : list
        List with the name of replicates of each condition (t_1_r_1, t_1_r_2,
        t_2_r_1, t_2_r_2)
    type_data : str
        Type of data. Type_data can be "counts" or "RPM".
    type_tables:
        Type of tables to join. type_tables can be "replicates" or "conditions".
    mode : str
        Method by which replicates of the same condition will join if
        type_tables is "replicates". If mode = "inner" the INNER JOIN method
        will be used, while if mode = "outer" the FULL OUTER JOIN method will
        be used. This parameter defaults to "inner".
    '''


    # Variables
    query_list = []
    sample_dict = {}
    if mode == 'outer':
        mode = ' FULL OUTER JOIN '
    elif mode == 'inner':
        mode = ' INNER JOIN '

    # Create dictionary with conditions names as keys and replicates list of
    # the condition as values.
    for sample in data_in:
        # Create condition's name
        sample_elements = sample.split("_")
        del sample_elements[3]
        new_sample_name = '_'.join(sample_elements)
        # Count the number of replicates of each condition
        if new_sample_name in sample_dict:
            sample_dict[new_sample_name].append(sample)
        else:
            sample_dict[new_sample_name] = [sample]   

        
    ### 1. MERGE REPLICATE TABLES
    ###########################################################################
    if type_tables == 'replicates':

        for key in sample_dict:

            # Samples list of the same condition
            sample_list = sample_dict[key]
            sample_list.sort()

            # New table = condition name
            new_table = key

            # Create query
            for j in range(len(sample_list)):
                if j == 0:
                    # First iteration
                    create = f'CREATE TABLE {new_table} AS '
                    select_seq = 'SELECT COALESCE(t1.seq' 
                    select_counts = f't1.{type_data} AS {sample_list[j]}'
                    on_section = '(t1.seq'
                    from_q = f' FROM {sample_list[j]} t1'
                    join = str()
                else:
                    # Complete query sections with the info from the rest of the replicates
                    select_counts = f'{select_counts}, t{str(j + 1)}.{type_data} AS {sample_list[j]}'
                    select_seq = f'{select_seq}, t{str(j + 1)}.seq'
                    join = f'{join}{mode}{sample_list[j]} t{str(j + 1)} ON {on_section}) = t{str(j + 1)}.seq'
                    # If second iteration
                    if j == 1:
                        on_section = f'COALESCE(t1.seq, t{str(j + 1)}.seq'
                    else:
                        on_section = f'{on_section}, t{str(j + 1)}.seq'
                   
        
            # Save query in query_list
            query = f'{create}{select_seq}) AS seq, {select_counts}{from_q}{join}'
            query_list.append(query.strip() + ";")

            # Remove tables used for joining
            for old_table in sample_list:
                drop_query = f'DROP TABLE {old_table};'
                query_list.append(drop_query)

            # Create Index of the table
            query_list.append(f'CREATE INDEX idx_{new_table} ON {new_table} (seq);')

    ### 2. MERGE CONDITION TABLES
    ###########################################################################      
    else:
        # First step in the query build
        first = True
        first_on = True

        # Indicates the current condition
        count_key = 1

        # Iterate conditions-replicates dictionary (t_1_r: t_1_r_1, t_1_r_2)
        for key in sample_dict:

            # Samples list of the same condition
            sample_list = sample_dict[key]
            sample_list.sort()

            # New table = condition name
            new_table = key

            # Final table
            new_table = "project_table"

            # Create query
            for i in range(len(sample_list)):

                # Write only when a new dictionary key is used and it is not the first one.
                if i == 0 and not first:
                    select_counts = f'{select_counts}, t{str(count_key + 1)}.{sample_list[i]} AS {sample_list[i]}'
                    select_seq = f'{select_seq}, t{str(count_key + 1)}.seq'
                    # First appearance of FULL OUTER JOIN in the query
                    if first_on:
                        join = f'{join} FULL OUTER JOIN {key} t{str(count_key + 1)} ON {on_section}) = t{str(count_key + 1)}.seq'
                        first_on = False
                    # Others
                    else:
                        join = f'{join} FULL OUTER JOIN {key} t{str(count_key + 1)} ON COALESCE{on_section}) = t{str(count_key + 1)}.seq'
                    on_section = f'{on_section}, t{str(count_key + 1)}.seq'
                    count_key += 1
                    
                # If it is the first step in the query build
                elif first:
                    create = f'CREATE TABLE {new_table} AS '
                    select_seq = 'SELECT COALESCE(t1.seq' 
                    select_counts = f't1.{sample_list[i]} AS {sample_list[i]}'
                    on_section = '(t1.seq'
                    from_q = f' FROM {key} t1 '
                    join = str()
                    first = False
                    
                # All other cases
                else:
                    select_counts = f'{select_counts}, t{str(count_key)}.{sample_list[i]} AS {sample_list[i]}'
        
        # Save query in query_list
        query = f'{create}{select_seq}) AS seq, {select_counts}{from_q}{join}'
        query_list.append(query.strip() + ";")

         # Create Index of the table
        query_list.append(f'CREATE INDEX idx_{new_table} ON {new_table} (seq);')


    ### 3. EXECUTE SQLITE QUERIES
    ###########################################################################

    # Connect to database
    sqliteConnection, cursor = connect_to_database(database_name)
    
    try:
        for query in query_list:
            cursor.execute(query)
        exit = False

    except Error as error:
        print("ERROR: Something went wrong during the join of the tables.")
        print(f'ERROR:\n{error}\n')
        exit = True
    
    finally:
        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()

        # If it fails, exit the program
        if exit:
            sys.exit()


def rep_counts_avg(database_name: str, table: str, new_table: str) -> str:
    """
    This function calculates the average counts (avg) of each sequence 
    belonging to the same condition (e.g. Control or Treated), that is,
    the average of the replicates (e.g. control1, control2). To do this, it
    accesses a table contained in a specified database, calculates the average
    of the counts of each sequence, and generates a table in the same database
    with the result (table_average).

    Parameters
    ----------
    database_name : str
        Absolute database path
    table : str
        Name of the input table.

    Returns
    -------
    new_table : str
        Name of the table containing the average of the counts of each sequence
        in the different conditions.
    """
        
    ## Connect to database
    sqliteConnection, cursor = connect_to_database(database_name)

    ## Get the counts table
    try:
        # Execute query
        table_sql = cursor.execute(f'SELECT * FROM {table}')

    except Error as error:
        # Error
        print(f'ERROR:\n{error}\n')

        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()

        # If it fails, exit the program
        sys.exit()
    
    ## Get condition replicates columns
    dic_colnames = {}
    seq_index = []
    for i in range(len(table_sql.description)):
        # Get column name
        col = table_sql.description[i][0]
        # Discard sequence columns
        if col[:3] == "seq":
            seq_index.append(i)
            continue
        # Remove replicate from the name
        col_elements = col.split("_")
        del col_elements[3]
        new_col = "_".join(col_elements)

        # Store in dictionary the columns of each condition (column = replicate)
        if new_col not in dic_colnames:
            dic_colnames[new_col] = [col]
        else:
            dic_colnames[new_col].append(col)

    ## Build Query
    final_columns_query = str()
    first_con = True

    # Iterate condition dictionary
    for con in dic_colnames:
        count_reps = int()
        new_colum_query = str()

        # Iterate replicates of each condition
        for j in range(len(dic_colnames[con])):
            count_reps += 1
            if j == 0:
                new_colum_query += f'CAST(({dic_colnames[con][j]}'
            else:
                new_colum_query += f' +  {dic_colnames[con][j]}'

        # new_column_query example: (control1 + control2) / 2 AS control
        new_colum_query += f') AS REAL ) / {str(count_reps)} AS {con}'

        if first_con:
            final_columns_query += new_colum_query
            first_con = False
        else:
            final_columns_query += ', ' + new_colum_query

    # Create final query
    query = f'CREATE TABLE IF NOT EXISTS {new_table} AS SELECT seq, {final_columns_query} FROM {table}'

    ## Execute final query
    try:
        # Execute query
        cursor.execute(query)
        exit = False

    except Error as error:
        print(f'ERROR:\n{error}\n')
        exit = True

    finally:
        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()
        
        # If it fails, exit the program
        if exit:
            sys.exit()

    return new_table


### 2. FORMAT CONVERSION FUNCTIONS

def sql_to_csv(database: str, table: str, path_out: str, columns: list,
               red_sample_dic: dict) -> None:
    """
    This function writes counts tables from a specific Sqlite database in .csv
    file.

    Parameters
    ----------
    database : str
        Absolute database path
    table : str
        Database table to write in .csv file.
    path_out : str
        Absolute path of the output file.
    columns : list
        List of table columns to be written to the new .csv file.
    red_sample_dic:
        Dictionary containing the original sample names and their associated
        reduced names.
    """

    ## 1. Connect to database
    sqliteConnection, cursor = connect_to_database(database)

    ## 2. Create query
    # Create a string with the columns to write from the table
    columns_str = 'seq'
    for col in columns:
        columns_str += f', ifnull({col}, 0) AS {col}'

    # Complete query with columns_str
    query = f'SELECT {columns_str} FROM {table}'

    ## 3. Execute query
    try:
        cursor.execute(query)
        exit = False

    except Error as error:
        print('ERROR. Something went wrong with the creation of the new .csv file')
        print(f'ERROR:\n{error}\n')
        exit = True  

    else:
        ## 4. Write table in .csv file
        print('Exporting data into CSV...')
        sys.stdout.flush()
        with open(path_out, 'w') as csv_file:
            csv_writer = csv.writer(csv_file)

            # Get original columns names
            original_columns = ['seq']
            for col in columns:
                red_name = col # p.e. t_1_r_3
                sample_name = red_sample_dic[red_name] # p.e. control_salt_1_5d
                original_columns.append(sample_name)

            # Write columns names and data
            csv_writer.writerow(original_columns)
            csv_writer.writerows(cursor)
    
    finally:    
        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()
        
        # If it fails, exit the program
        if exit:
            sys.exit()
    

def fasta_to_csv(fasta_file: str, sep: str, new_file_path: str) -> None:
    '''
    This function creates a file with two columns (sequence, identifier)
    delimited by a specific separator (sep) from a .fasta file. 

    Parameters
    ----------
    fasta_file : str
        Absolute path of the .fasta file.
    sep : str
        Separator delimiting the sequence and identifier columns in
        the new file.
    new_file_path : str
        Absolute path of the new file.
    '''
    # Parse .fasta file
    fasta_generator = parse_fasta_file(fasta_file)
    
    # Write in the new file
    with open(new_file_path, 'w') as tsv_file:
        tsv_file.write(f'seq{sep}header\n') 
        for line in fasta_generator:
            try: 
                header, seq = line
                tsv_file.write(f'{seq}{sep}{header}\n')
            except StopIteration as error:
                print('Stop Iteration occured')
                print(error)
                break


### 3. FILTERING FUNCTIONS

def filter_by_rnacentral(library_list: list, Rnacentral_dir: str,
                        path_write: str) -> None:
    '''
    This function filters the libraries by eliminating those sequences that
    map against the RNAcentral database.

    Parameters
    ----------
    library_list : list
        List with the absolute paths of the libraries to be filtered
        with rnacentral.
    RNAcentral_database : str
        Path of the directory where the RNAcentral database is located.
    path_write : str
        Path of the directory where all the libraries filtered by RNAcentral
        will be stored.

    '''

    # Create the directory where the filtered libraries will be stored
    path_align = f'{path_write}/alignments'
    os.system(f'mkdir -p {path_align}')
    
    # Save current directory
    script_path = os.getcwd()

    # Assign this directory as the current directory
    os.chdir(path_align)

    # Create directory where the indexed database will be stored
    path_index = f'{path_align}/Index'
    os.system(f'mkdir -p {path_index}')
    
    print('\nPercentage of sequences aligned on RNAcentral database:')

    # Iterate library list
    for library in library_list:
        
        # Get the file name
        name_filtered = os.path.basename(library) # SRRXXXXXXX_si_filtered.fasta
    
        # Split the name and select the SRR
        name = name_filtered.split('_')[0] # SRRXXXXXXX

        # Execute Bowtie
        os.system(f'bowtie -x {Rnacentral_dir}/Index/rnacentral --best -v 0 -k 1 --norc -f {library} -S {name}.sam --al {name}_RNAcentral_aligned.fasta --un ./../{name}_RNAcentral_filtered.fasta')
    
        # Number of reads that have aligned with the database
        try:
            # Open file
            align = open(f'{path_align}/{name}_RNAcentral_aligned.fasta', 'r')
            
        except FileNotFoundError:
            # If no reads have aligned, the file will not exist.
            n1 = 0
            pass
        else:
            # Count the number of sequences
            n1 = len([line for line in align if line.strip()[0] == '>'])
        finally:
            align.close()
        
        # Number of reads that have NOT aligned with the database
        try:
           unalign = open(f'{path_write}/{name}_RNAcentral_filtered.fasta', 'r')
        
        except FileNotFoundError:
            # If all reads have aligned, the file will not exist.
            n2 = 0
            pass
        else:
            # Count the number of sequences
            n2 = len([line for line in unalign if line.strip()[0] == '>'])
        finally:
            unalign.close()

        # Get the percentage of sequences that have aligned with the database.
        try:
            total = n1 + n2
            percentage_aligned = (n1/total)*100
        except ZeroDivisionError:
            # Total should not be 0. This only occurs if the files have not
            # been generated, so they cannot be opened.
            print('ERROR: there has probably been a problem with the bowtie alignment.')
            sys.exit()
        else:
            print(f'--{name}: {str(percentage_aligned)}%')
        
    # If no sequence has aligned with RNAcentral, delete the alignment directory.
    if percentage_aligned == 0:
        os.system('rm -R ../alignments')

    # If alignment has occurred, delete everything except the files containing
    # the aligned sequences.
    else:
        os.system('pwd')
        os.system('rm -r ./Index/')
        os.system('rm ./*.sam')

    os.chdir(script_path)


### 4. PARALLELIZED FUNCTIONS

def all_absolute_counts_processing (libraries_list: list, path_read: str,
                                 path_write: str, num_process: int) -> None:
    '''
    This function is used to parallelize the all_absolute_counts function,
    which calculates the absolute counts of each library

    Parameters
    ----------
    libraries_list : list
        List of input fasta file names (libraries). As it is parallelized,
        each process will call this function with a different list of fasta
        files. 
    path_read : str
        Absolute path of the directory where the input libraries are located.
    path_write : str
        Absolute path of the directory where the absolute count tables will
        be stored.
    num_process : int
        Number of processes
    '''
    number_files = len(libraries_list)

    # If the number of libraries is <= to the maximum number of processes (30)...
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
        
        if i == num_process-1:
            pf = number_files
        else:
            pf = pi + distribution
        
        processes.append(multiprocessing.Process(target = all_absolute_counts ,
                                                args=(libraries_list[pi:pf],
                                                      path_read, path_write))) 
        processes[i].start()
        print(f'Process {i} launched.')
        
        pi = pf
        
    for p in processes:
        p.join()
    
    print('\nDONE!')
    
    
def all_absolute_counts (libraries_list: list, path_read: str,
                         path_write: str) -> None:
    '''
    This function calculates the absolute counts of a libraries list.

    Parametros
    ----------
    libraries_list : list
        List of input fasta file names (libraries). As it is parallelized,
        each process will call this function with a different list of fasta
        files. 
    path_read : str
        Absolute path of the directory where the input libraries are located.
    path_write : str
        Absolute path of the directory where the absolute count tables will
        be stored.
    '''
    os.chdir(path_write)

    # Iterate the library fasta files list.
    for name in libraries_list:

        # Create the name of the new file
        if name.endswith('.fa'):
            run_name = name.replace('.fa', '')
        elif name.endswith('.fasta'):
            run_name = name.replace('.fasta', '')

        abs_name = f'{run_name}_abs.csv'  # SRRXXXXXXX_abs_counts.csv

        # Absolute input and output paths
        path_file_in = f'{path_read}/{name}'
        path_file_out = f'{path_write}/{abs_name}'
        
        # Create the count table using bash. It is faster
        command = f'echo \"seq,counts\" > {path_file_out}'
        os.system(command)
        command = "grep -v '>' " + path_file_in + " | sort | uniq -c | sort -nr | awk 'BEGIN{FS=\" \"; OFS=\",\"} {print $2, $1}' >> " + path_file_out
        os.system(command)


def all_rpm_counts_processing(abs_list: list, path_read: str, path_write: str,
                              num_process: int) -> None:
    '''
    This function is used to parallelize the all_rpm_counts function, which
    calculates the reads per millon (RPM) of each library.

    Parameters
    ----------
    abs_list : list
        List of input csv file names. These files contain the absolute count
        tables. As it is parallelized, each process will call this function
        with a different list of csv files. 
    path_read : str
        Absolute path of the directory where the input csv files are located.
    path_write : str
        Absolute path of the directory where the RPM tables will be stored.
    num_process : int
        Number of processes.
    '''

    number_files = len(abs_list)

    # If the number of libraries is <= to the max number of processes (30)...
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
        
        if i == num_process-1:
            pf = number_files
        else:
            pf = pi + distribution
        
        processes.append(multiprocessing.Process(target = all_rpm_counts ,
                                                 args=(abs_list,
                                                       path_read,
                                                       path_write) )) 
        processes[i].start()
        print(f'Process {i} launched.')
        
        pi = pf
        
    for p in processes:
        p.join()
    
    print('\nDONE!')


def all_rpm_counts (abs_list: list, path_read: str, path_write: str) -> None:
    '''
    This function calculates the reads per millon (RPM) of a libraries list.

    Parameters
    ----------
    abs_list : List
        List of input csv file names. These files contain the absolute count
        tables. As it is parallelized, each process will call this function
        with a different list of csv files. 
    path_read : String
        Absolute path of the directory where the input csv files are located.
    path_write : String
        Absolute path of the directory where the RPM tables will be stored.
    '''

    # Iterate the absolute counts files list.
    for abs_name in abs_list:
        
        # Create the name of the new file
        run_name = abs_name.replace('_abs.csv', '')
        rpm_name = f'{run_name}_rpm.csv'

        # Absolute input and output paths
        path_file_out = f'{path_write}/{rpm_name}'
        path_file_in = f'{path_read}/{abs_name}'

        # Calculate the total number of reads
        chunksize = 1000
        sum_reads = 0
        with open(path_file_in) as reader_sum:
            sum_df = pd.read_csv(reader_sum, sep=",", chunksize=chunksize, low_memory=False)
            for chunk in sum_df:
                sum_reads += chunk['counts'].sum()
        
        # Create the RPM table
        with open(path_file_in, "r") as reader, open(path_file_out, 'w') as writer:
            count_tab = pd.read_csv(reader, sep=",", chunksize=chunksize, low_memory=False)
            first_chunk = True
            for chunk in count_tab:
                # Insert chunk in a table and calculate RPM
                chunk_df = pd.DataFrame(chunk, columns=['seq', 'counts'])
                abs_ = np.array(list(chunk_df['counts']), dtype = float) 
                rpm = (abs_ * 1000000) / sum_reads
                chunk_df.insert(value = rpm, column='RPM', loc=2)
                chunk_df = chunk_df.iloc[:,[0,2]]
                chunk_df.columns = ['seq', 'RPM']
                # Save result in a new file
                if first_chunk:
                    chunk_df.to_csv(writer, index = None, mode='a') 
                    first_chunk = False
                else:
                    chunk_df.to_csv(writer, index = None, mode='a', header = False)


### 5. METADATA FUNCTIONS

def get_sample_metadata(path: str, sample: str) -> list:
    '''
    This function extracts the sample metadata from the metadata table of the
    project to which it belongs. 

    Parameters
    ----------
    path : str
        Absolute path of the txt file containing the project metadata.
    sample : str
        Sample name (e.g. SRRXXXX)

    Returns
    -------
    SRR_metadata : list
        Metadata list of the sample of interest.
    '''
   
    SRR_metadata = []
    try:
        # Read metadata csv file
        metadata_table = genfromtxt(path, delimiter=',',
                                    dtype=None, encoding=None)
    except Exception as e:
        print(f'Unable to read metadata table of {sample}')
        print(f'Exception: {e}')
        sys.exit()
    else:
        # Iterate the metadata table
        for row in range(len(metadata_table)):
            # Get SRR and remove spaces at beginning and end of string
            run = metadata_table[row][2].strip()
            # Find the row in which the sample is located
            if run == sample:
                # If the SRRs match, save their metadata in a list.
                SRR_metadata = metadata_table[row]
                break
    
    return SRR_metadata


def get_informative_variables(path_metadata: str) -> list:
    '''
    This function identifies informative variables within a metadata table,
    understanding informative variables as those with more than one level
    (e.g. Time: 0h, 3h and 6h). 

    Parameters
    ----------
    path_metadata : str
        Absolute path of the csv containing the metadata of a specific project.

    Returns
    -------
    inf_var_index : list
        List with the indexes or positions (columns) in which the information
        of the informative variables is located (p.e. stress -> 4).
    '''
    try:
        # Get metadata
        metadata_table = genfromtxt(path_metadata, delimiter=',',
                                    dtype=None, encoding=None,)
    except Exception as e:
        print(f'Unable to obtain informative variales from: \
              {os.path.basename(path_metadata)}')
        print(f'Exception: {e}')
        sys.exit()
    else:
        # Create dataframe
        df = pd.DataFrame(metadata_table)

        # List of columns or indexes in which the informative variables are
        # found within the data.frame. Indexes 3, 4 and 5 are condition
        # (cntl vs trtd), stress (p.e. cold) and replica (p.e. 1,2 or 3),
        # which will always be considered.
        inf_var_pos = [3, 4, 5]
        
        # Get the number of columns except the first 6 (already saved or
        # not important).
        ncol = len(df.columns) - 6
        
        # Iterate ncol
        for col in range(ncol):
            # Extract the "col" column
            col_list =  list(df.iloc[:,col + 6])

            # Check if the variable contains different levels
            niveles = pd.unique(col_list[1:])

            # Save the indexes of those variables that have more than 1 level.
            if len(niveles) != 1 :
                indice = col + 6
                inf_var_pos.append(indice)
    
        return inf_var_pos


### 6. OTHER FUNCTIONS

def parse_fasta_file(fasta_file: str) -> None:
    '''
    This function parse a fasta file returning a generator object containing
    the headers and the sequences. If this object is iterated outside the
    function, we can access a new sequence and its header at each iteration
    of the loop.

    Parameters
    ----------
    fasta_file : str
        Absolute path of the fasta file
    '''
    try:
        with open(fasta_file) as fasta:
            iterator = (x[1] for x in itertools.groupby (fasta, lambda line: line[0] == '>'))
            for header in iterator:
                # Drop the ">"
                headerStr = header.__next__()[1:].strip()
                # join all sequence lines to one
                seq = ''.join(s.strip() for s in iterator.__next__())
                yield(headerStr,seq)

    except StopIteration:
        print('Error. Check that the files entered are fasta')
        sys.exit()


def fusion_tables (list_abs: list, list_rpm: list, path_write: str,
                  path_metadata: str, mode: str) -> None:
    '''
    This function merges all absolute count tables into a single table and all
    RPM tables into another different table. In addition, it also creates a
    table for both types of counts (Abs and RPM) by calculating the average
    of the replicates.

    Parameters
    ----------
    list_abs : list
        Absolute paths list of the absolute count table for each library.
    list_rpm : str
        Absolute paths list of the RPM table for each library.
    path_write : str
        Path to the directory where the tables generated in this function
        will be stored.
    path_metadata : str
        Absolute project metadata path
    mode : str
        Method used for the replicates joining. It can be "inner" (PCA) or
        "outer" (differential expression analysis).
    '''

    # Variables
    red_sample_list = []
    red_sample_dic = {}
    red_condition_dic = {}
    con_num_dic = {}
    con_count = 1

    # Get project name
    project_name = list_abs[0].split("/")[-2].split('_')[0]

    # Create temporary directory
    os.system(f'mkdir -p ./tmp_{project_name}')

    # Get informative variables from the project metadata (columns from 0-n)
    inf_var = get_informative_variables(path_metadata)

    # Iterate a number of times equal to len(list_abs) and len(list_rpm)
    for index1 in range(len(list_abs)):

        ## 1. BUILD SAMPLE NAME (e.g. control_cold_2_Col-0)
        #######################################################################

        # If it has been filtered by RNAcentral...
        if list_abs[index1].endswith('_RNAcentral_filtered_abs.csv'):
            
            # Get project and run name from the path
            path_elements = list_abs[index1].split('/')
            run_name = path_elements[-1].replace('_RNAcentral_filtered_abs.csv', '')
        
        # If it has NOT been filtered by RNAcentral
        else:
            # Get project and run name from the path
            path_elements = list_abs[index1].split('/')
            run_name = path_elements[-1].replace('_tr_abs.csv', '')
        
        print(f'Obtaining metadata of sample: {run_name}', end='')
        sys.stdout.flush()

        # Get SRR metadata
        metadata = get_sample_metadata(path_metadata, run_name)
        
        # If there is no sample metadata, exit.
        if len(metadata) == 0:
            print(f'\nERROR:  {run_name} run does not have metadata.\n')
            return
      
        # Build the sample name with the metadata using the informative variables
        sample_name = str()
        con_name = str()
        first_var = True
        for index2 in inf_var:
            # Get the name of the informative variable (e.g. control or treated)
            var = str(metadata[index2]).strip()

            # Add it to the final name of the sample
            if first_var:
                sample_name += var
                con_name += var
                first_var = False

            else:
                # If it is not the first variable, we precede it with a "_".
                sample_name += '_' + var # e.g. control_drought_1_5d
                if index2 != 5:
                    con_name += '_' + var # e.g. control_drought_5d
                if index2 == 5:
                    replicate_num = var # e.g. replicate = 1
        
        # Save in dictionary the condition and its identifier number.
        if con_name not in con_num_dic:
            con_num_dic[con_name] = con_count # p.e. control = 1, treated = 2
            condition_num = con_count
            con_count += 1
        else:
            condition_num = con_num_dic[con_name]
        
        # Generate shortened name and establish correspondence with original name
        # Replicate name. Necessary to create the abs and RPM tables.
        new_red_name = f't_{str(condition_num)}_r_{str(replicate_num)}' # p.e t_1_r_3
        red_sample_dic[new_red_name] =  sample_name # Replicate
        red_sample_list.append(new_red_name)
        print(f' ({new_red_name})')

        # Condition name. Necessary to create the tables with the abs and RPM averages
        new_red_con = f't_{str(condition_num)}_r'
        if new_red_con not in red_condition_dic:
            red_condition_dic[new_red_con] =  con_name
        
        ## 2. LOAD ABSOLUTE COUNTS TABLE INTO CORRESPONDING SQLITE DATABASE
        #######################################################################
        # Extract name of the database to be created (project)
        db_abs_name = f'./tmp_{project_name}/{list_abs[index1].split("/")[-2]}_{mode}.db' # ej. PRJNA277424_abs_RF.db
        # Insert absolute counts of each sample in db (1 sample = 1 table in db)
        print(f'Inserting {sample_name} table in database (Absolute Counts)...')
        sys.stdout.flush()
        insert_to_database(db_abs_name, new_red_name, list_abs[index1], 'counts_abs')

        ## 3. LOAD RPM TABLE INTO CORRESPONDING SQLITE DATABASE
        #######################################################################
        # Extract name of the database to be created (project)
        db_rpm_name = f'./tmp_{project_name}/{list_rpm[index1].split("/")[-2]}_{mode}.db' # ej. PRJNA277424_rpm_RF.db
        # Insert RPM of each sample in db (1 sample = 1 table in db)
        print(f'Inserting {sample_name} table in database (RPM)...')
        sys.stdout.flush()
        insert_to_database(db_rpm_name, new_red_name, list_rpm[index1], 'counts_rpm')
        
    
    ## 4. JOIN THE DIFFERENT REPLICATES OF EACH CONDITION IN A SINGLE TABLE
    ###########################################################################
    print('Joining samples of the same condition (Replicates)...')
    sys.stdout.flush()

    # Tuple with the arguments of the merge_counts_tables function for Abs and RPM.
    list_args_merge_rep = [(db_abs_name, red_sample_list, 'counts', 'replicates', mode),
                        (db_rpm_name, red_sample_list, 'RPM', 'replicates', mode)]          
    merge_counts_tables_processing(list_args_merge_rep)
    print('DONE!\n')
    sys.stdout.flush()
 

    ## 5. JOIN ALL THE PROJECT COUNTS IN A SINGLE TABLE
    ###########################################################################

    print('Joining condition tables...')
    sys.stdout.flush()

    # Tuple with the arguments of the merge_counts_tables function for Abs and RPM.
    list_args_merge_con = [(db_abs_name, red_sample_list, 'counts', 'conditions', 'outer'),
                        (db_rpm_name, red_sample_list, 'RPM', 'conditions', 'outer')]
                    
    merge_counts_tables_processing(list_args_merge_con)
    
    print('DONE!\n')
    sys.stdout.flush()
    

    ## 6. WRITE THE RESULTING TABLE IN A .CSV FILE
    ###########################################################################

    # Output paths
    path_write_abs = f'{path_write}/fusion_abs-{mode}.csv'
    path_write_rpm = f'{path_write}/fusion_rpm-{mode}.csv'

    # Sort sample list
    red_sample_list.sort()

    # Abs counts
    final_table = 'project_table'
    print('Saving tables resulting from the joining process in a csv file (Absolute Counts)...')
    sql_to_csv(db_abs_name, final_table, path_write_abs, red_sample_list, red_sample_dic)

    # RPM
    print('Saving tables resulting from the joining process in a csv file (RPM)...')
    sql_to_csv(db_rpm_name, final_table, path_write_rpm, red_sample_list, red_sample_dic)
    sys.stdout.flush()


    ## 7. CALCULATE REPLICATES AVERAGE AND SAVE IT IN .CSV FILE
    ###########################################################################

    # Out paths
    path_write_mean_abs = f'{path_write}/fusion_abs_mean-{mode}.csv'
    path_write_mean_rpm = f'{path_write}/fusion_rpm_mean-{mode}.csv'

    ## Calculate average
    # Abs counts
    print('Calculating average of replicates of each condition (Absolute Counts)...')
    sys.stdout.flush()
    rep_counts_avg(db_abs_name, final_table, final_table + '_avg')

    # RPM
    print('Calculating average of replicates of each condition(RPM)...')
    sys.stdout.flush()
    rep_counts_avg(db_rpm_name, final_table, final_table + '_avg')
    print("DONE!\n")
    sys.stdout.flush()

    ## Save in .csv
    list_condition = list(red_condition_dic.keys())
    # Abs counts
    print('\nSaving AVG-Abs tables in a csv file...')
    sys.stdout.flush()
    sql_to_csv(db_abs_name, final_table + '_avg', path_write_mean_abs, list_condition, red_condition_dic)

    # RPM
    print('\nSaving AVG-RPM tables in a csv file...')
    sys.stdout.flush()          
    sql_to_csv(db_rpm_name, final_table + '_avg', path_write_mean_rpm, list_condition, red_condition_dic)
    print('Files saved!\n')
    sys.stdout.flush()

    ## 8. DELETE DATABASES AFTER USE
    ###########################################################################
    os.system(f'rm -r ./tmp_{project_name}')

    
## MAIN PROGRAM

def main():
    '''
    Main program
    '''
    parser = argparse.ArgumentParser(prog='FILTAB V2', 
                                     description='''This program filters sRNA \
                                        sequences using RNAcentral to discard \
                                        rRNA, tRNA, snRNA and snoRNA sequences \
                                        present in the libraries of a specific \
                                        project, generates tables of absolute \
                                        counts and reads per million for each \
                                        of them, and then joins these tables \
                                        to form a single table for the project.''',  
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-p', '--path-project', type=str, nargs=1,  
                        help='Absolute path of project libraries folder.')
    parser.add_argument('-r', '--path-results', type=str, nargs=1, 
                        help='Absolute path of results folder.')
    parser.add_argument('-t', '--threads', type=int, nargs='?', const=1, default=1,
                        help='Number of threads (default: %(default)s)')
    parser.add_argument('-x', '--path-metadata', type=str, nargs=1,  
                        help='Absolute path of metadata samples table.')
    parser.add_argument('-z', '--path-rnacentral', type=str, nargs=1,  
                        help='Absolute path of RNAcentral database.')

    parser.add_argument('--version', action='version', version='%(prog)s 2.0')
    
    args = parser.parse_args()
    
    
    ###########################################################################
    #                        1. CHECK ARGUMENTS                               #
    ###########################################################################
    
    ## 1.1. CHECK ARGUMENTS
    ###########################################################################

    try:
        path_project_dir = args.path_project[0]
        path_res = args.path_results[0]
        num_process = args.threads
        path_metadata = args.path_metadata[0]
    except:
        print('ERROR: You have inserted a wrong parameter or you are missing a parameter.')
        parser.print_help()
        sys.exit()
    

    ## 1.2. CHECK RNACENTRAL ARGUMENT
    ###########################################################################
    try:
        path_RNAcentral = args.path_rnacentral[0]
        R = 1
    except:
        print('\nIMPORTANT: No RNAcentral database has been inserted.')
        R = 0
    else:
        print('\nIMPORTANT: RNAcentral database has been inserted.')

        # Get the database directory path
        RNAcentral_dir = os.path.dirname(path_RNAcentral)
    
        # Check if RNAcentral has already been indexed.
        rnacentral_list = os.listdir(RNAcentral_dir)
        if 'Index' in rnacentral_list:
            print('RNAcentral was previously Indexed!')
        else:
            # Indexar RNAcentral
            print('\nRNAcentral Database is not indexed.\n')
            print('INDEXING RNAcentral DATABASE...\n')
            path_index = f'{RNAcentral_dir}/Index/'
            os.system(f'mkdir -p {path_index}')
            os.system(f'bowtie-build --threads ' + str(num_process) +
                      ' --large-index ' + path_RNAcentral + ' ' + path_index +
                      'rnacentral')

    #######################################################################
    #                   4. CALCULATE ABSOLUTE COUNTS AND RPM              #
    #######################################################################
    
    # Get the path of the directory in which the trimmed libraries are stored
    path_project_split = path_project_dir.split("Libraries")
    path_lib = path_project_split[0] + "Libraries"

    # Get the project and the species name
    project = os.path.basename(path_project_dir)
    species = os.path.basename(os.path.dirname(path_project_dir))
    
    if R == 1:
        
        ### 4.A.1. FILTER BY RNACENTRAL
        ###################################################################

        ## Input path
        path_read = path_project_dir

        ## Set suffix for directories
        suffix = '_RF'

        ## Output path
        folder_write = '03-RNAcentral_filtered'
        path_write = f'{path_lib}/{folder_write}/{species}{suffix}/{project}{suffix}'
    
        ## List of input fasta files
        fasta_files_list = []
        fasta_files = os.listdir(path_read)
        for file in fasta_files:
            path_fasta_files = f'{path_read}/{file}'
            if file.endswith('.fasta') or file.endswith('.fa'):
                fasta_files_list.append(path_fasta_files)
        
        ## Filter by RNAcentral
        print('\nFILTER BY RNAcentral DATABASE.\n')
        print('Filtering...\n')
        sys.stdout.flush()
        filter_by_rnacentral(fasta_files_list, RNAcentral_dir, path_write)
    
    elif R == 0:

        ### 4.A.2. DO NOT FILTER BY RNACENTRAL
        ###################################################################

        ## Set suffix for directories
        suffix = ''

        ## Input path for the following section
        path_write = path_project_dir


    ### 4.B. CALCULATE ABSOLUTE COUNTS
    #######################################################################
    
    ## Create new output path (absolute counts)
    folder_write_abs = f'01-Absolute_counts{suffix}'
    path_write_abs = f'{path_res}/{folder_write_abs}/{species}_abs{suffix}/{project}_abs{suffix}'

    ## Get input path
    path_read_abs = path_write

    ## List of input fasta files
    fasta_read_abs_list = os.listdir(path_read_abs)
    for file in fasta_read_abs_list:
        if file.endswith('.fasta') or file.endswith('.fa'):
            continue
        else:
            fasta_read_abs_list.remove(file)
    
    ## Create output directory
    print('\nCREATING THE ABSOLUTE COUNT TABLE OF EACH LIBRARY.\n')
    print('Building the directory...\n')
    sys.stdout.flush()
    try:
        os.makedirs(path_write_abs)
    except FileExistsError:
        print(f'The folder {path_write_abs.split("/")[-2]} already exists.')
        pass
    
    ## Get absolute counts for each library
    print('\nCalculating the absolute counts...')
    sys.stdout.flush()
    t_ini = time()
    all_absolute_counts_processing(fasta_read_abs_list, path_read_abs, path_write_abs, num_process)
    t_end = time()
    print('Time: %.3f s' % (t_end - t_ini))
    sys.stdout.flush()

    
    ### 4.C. CALCULATE RPM
    #######################################################################
    
    ## Create new output path (RPM)
    folder_write_rpm = f'02-RPM_counts{suffix}'
    path_write_rpm = f'{path_res}/{folder_write_rpm}/{species}_rpm{suffix}/{project}_rpm{suffix}'

    ## Get input path (absolute counts)
    path_read_rpm = path_write_abs

    ## List of absolute counts files (input)
    abs_csv_list = os.listdir(path_read_rpm)

    ## Create output directory
    print('\nCREATING THE RPM TABLE OF EACH LIBRARY.\n')
    print('Building the directory...\n')
    sys.stdout.flush()
    
    try:
        os.makedirs(path_write_rpm)
    except FileExistsError:
        print(f'The folder {path_write_rpm.split("/")[-2]} already exists.')
        pass
    
    ## Get RPM for each library
    print('\nCalculating the RPM counts...')
    sys.stdout.flush()
    t_ini = time()
    all_rpm_counts_processing(abs_csv_list, path_read_rpm, path_write_rpm, num_process)
    t_end = time()
    print('Time: %.3f s' % (t_end - t_ini))
    sys.stdout.flush()
        
    
    ### 4.D. JOIN ABSOLUTE COUNTS AND RPMS TABLES OF THE DIFFERENT SAMPLES
    #######################################################################
    
    ## Paths
    path_read_join_abs = path_write_abs
    path_read_join_rpm = path_write_rpm
    folder_write_join = f'03-Fusion_count_tables{suffix}'
    path_write_join = f'{path_res}/{folder_write_join}/{species}/{project}'
    
    ## List of the input csv files paths (absolute counts)
    csv_files_list_abs = []
    csv_files_abs = os.listdir(path_read_join_abs)
    for abs_file in csv_files_abs:
        path_csv_abs_file = f'{path_read_join_abs}/{abs_file}'
        if abs_file.endswith('.csv'):
            csv_files_list_abs.append(path_csv_abs_file)
    
    ## List of the input csv files paths (RPM)
    csv_files_list_rpm = []
    csv_files_rpm = os.listdir(path_read_join_rpm)
    for rpm_file in csv_files_rpm:
        path_csv_rpm_file = f'{path_read_join_rpm}/{rpm_file}'
        if rpm_file.endswith('.csv'):
            csv_files_list_rpm.append(path_csv_rpm_file)
    
    ## Create output directory
    print('\nJOINING THE RPM AND ASBOLUTE TABLES IN TWO TABLES.\n')
    print('Building the directory...\n')
    sys.stdout.flush()
    try:
        os.makedirs(path_write_join)
    except FileExistsError:
        print(f'The folder {path_write_join.split("/")[-2]} already exists.')
        pass

    # Ruta de los metadatos del proyecto en cuestion
    path_metadata_abs = f'{path_metadata}/{species}_m_{project}.txt'

    # Join tables
    print('\nJoining the tables mode = outer between replicates...')
    sys.stdout.flush()
    fusion_tables (csv_files_list_abs, csv_files_list_rpm, path_write_join, path_metadata_abs, 'outer')

    ###############################################################
        
## CALL THE MAIN PROGRAM

if __name__ == '__main__':
    '''
    Call to the main program
    '''
    main()
