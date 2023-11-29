#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#******************************************************************************
#  
#   Divide_filter_by_subproject.py
#
#   This program divides the project count tables by subproject writing
#   the table section of each subproject in a separate file. To do this,
#   it uses the metadata of each project to group the columns of the
#   treatment samples with the columns of the control samples with which
#   they are associated. Finally, it extracts each group of columns from
#   the original table (previously inserted into an SQLite database) and
#   writes them to a new separate "CSV" file while selecting those
#   sequences that have 5 or more counts in at least 5 samples. In case
#   the subproject has less than 5 samples, all the subproject samples
#   must have at least 5 counts.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 21/09/2023
#   Version: 1.1
#
#
#******************************************************************************


## IMPORT MODULES

import argparse
import csv
import numpy as np
import os
import pandas as pd
import sqlite3
import sys


## FUNCTIONS

### 1. SQLITE FUNCTIONS

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
        
    except sqlite3.Error as error:
        print('ERROR. Unable to connect to the database')
        print(f'ERROR. {error}')
        sys.exit()

    else:
        return sqliteConnection, cursor


def insert_to_database (database_name: str, table_name: str, data_path: str, columns: list):
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
    columns : list
        Columns list of the table to be inserted
    '''

    # Connect to database
    sqliteConnection, cursor = connect_to_database(database_name)    

    # Create dictionary with original and temporal names of samples
    id_dic = {}
    cond_list = []
    count = 0
    for column in columns:
        if column == 'seq':
            continue
        # Get condition name and rep
        column_elements = column.split('_')
        rep = column_elements.pop(2)
        condition = '_'.join(column_elements)

        # Create sample id
        if condition not in cond_list:
            count += 1
            id_dic[column] = f't_{str(count)}_r_{rep}'
            cond_list.append(condition)    
        else:
            id_dic[column] = f't_{str(count)}_r_{rep}'
    
    # Create columns list and query section
    query_section = ''
    new_columns = []
    for i, column in enumerate(columns):
        if column == 'seq':
            query_section = '(seq TEXT, '
            new_columns.append(column)
        elif i == len(columns) - 1:
            query_section = f'{query_section}{id_dic[column]} INT);'
            new_columns.append(id_dic[column])
        else:
            query_section = f'{query_section}{id_dic[column]} INT, '
            new_columns.append(id_dic[column])
    
    # Read data to insert
    chunksize=1000000
    data_to_insert = pd.read_csv(data_path, sep=',',
                                 chunksize=chunksize,
                                 low_memory=False,
                                 header=None,
                                 names=new_columns)

    try:
        # Pragma adjust
        cursor.execute('PRAGMA synchronous = OFF')
        cursor.execute('PRAGMA journal_mode = OFF')

        # Create table
        cursor.execute(f'CREATE TABLE IF NOT EXISTS {table_name}{query_section}')

        # Insert chunks in table
        for chunk in data_to_insert:
            # Insert chunk data into table
            chunk.to_sql(name=table_name,
                         con=sqliteConnection,
                         if_exists='append',
                         index=False)       
        sys.stdout.flush()

        # Create table index to increase query speed
        sys.stdout.flush()
        cursor.execute(f'CREATE INDEX idx_{table_name} ON {table_name} (seq);')
        sys.stdout.flush()
        exit = False

    except sqlite3.Error as error:
        print('ERROR: Something went wrong during the data insertion.')
        print(error)
        exit = True
    
    finally:
        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()
        
        # If it fails, exit the program
        if exit:
            sys.exit()
        
        return id_dic
    
    
def sql_to_csv (database: str, table: str, path_out: str, columns: list, names_dic: dict) -> None:
    '''
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
    names_dic:
        Dictionary containing the original sample names and their associated
        reduced names.
    '''

    ## 1. Connect to database
    sqliteConnection, cursor = connect_to_database(database)

    ## 2. Create query
    # Create a string with the columns to write from the table
    columns_str = 'seq'
    for col in columns:
        columns_str += f', {names_dic[col]}'

    # Complete query with columns_str
    query = f'SELECT {columns_str} FROM {table}'

    ## 3. Execute query
    try:
        cursor.execute(query)
        exit = False

    except sqlite3.Error as error:
        print('ERROR. Something went wrong with the creation of the new .csv file')
        print('ERROR: ', end='')
        print(error)
        exit = True  

    else:
        ## 4. Write table in .csv file
        with open(path_out, 'w') as csv_file:
            csv_writer = csv.writer(csv_file)

            # Iterate SQLite cursor object          
            for i, line in enumerate(cursor):

                # Write header
                if i == 0:
                    csv_writer.writerow(line)
                    continue

                # Select sequences that have more than 5 counts in at least 5
                # sample
                count = len([i for i in line[1:] if i > 5])
                if count >= 5:
                    csv_writer.writerow(line)

                # If there are less than 5 samples, check that the condition is
                # fulfilled for all of them.
                elif count < 5 and count == len(line) - 1:
                    csv_writer.writerow(line)
    
    finally:    
        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()
        
        # If it fails, exit the program
        if exit:
            sys.exit()


### 2. METADATA FUNCTIONS
def get_project_metadata (path: str) -> pd.DataFrame:
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
        print(f'Exception: {e}')
        sys.exit()
    else:
        return df


def get_informative_variables (df: pd.DataFrame) -> list:
    '''   
    This function identifies informative variables within a metadata table,
    understanding informative variables as those with more than one level
    (e.g. Time: 0h, 3h and 6h).
    
    Parameters
    ----------
    df : _type_
        _description_

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
        col_list =  np.array(df.iloc[:,col + 6])

        # Check if the variable contains different levels
        niveles = pd.unique(col_list[1:])

        # Save the indexes of those variables that have more than 1 level.
        if len(niveles) != 1 :
            indice = col + 6
            inf_var_pos.append(indice)

    return inf_var_pos


### 3. DIVIDE BY SUBPROJECTS FUNCTION
def divide_project_by_subproject (project_samples, project_metadata):
    '''
    This function receives a list of samples from a project and groups them
    according to the subproject they belong to. It does this by comparing
    control samples with treatment samples to identify which treatment samples
    are associated with which control samples. Once identified, it groups them
    in a list together with their control samples. The function returns a list
    containing the list corresponding to each of these groups.
    
    Parameters
    ----------
    project_samples : list
        Project samples list
    project_metadata : str
        Absolute path of project metadata file

    Returns
    -------
    groups_list : list
        List of sample groups
    '''

    # Get project metadata
    meta_table = get_project_metadata(project_metadata)
    infvar = get_informative_variables(meta_table)

    # Get column names
    columns = [meta_table.columns[var] for var in infvar]
    del columns[2] # Delete replicate column

    # Required lists
    controls_condition = []
    treated_condition = []
    control_samples = []
    treated_samples = []

    # Separate controls and treatments into two separate lists
    for sample in project_samples:
        # Discard sequence column
        if sample == 'seq':
            continue
        # Delete replicate column
        sample_elements = sample.split('_')
        del sample_elements[2]
        # Save control samples
        if sample.startswith('control'):
            # Save control samples in list 
            control_samples.append(sample)
            # Save control conditions in list
            if sample_elements not in controls_condition:
                controls_condition.append(sample_elements)
                
        # Save treated conditions in list
        elif sample.startswith('treated'):
            treated_samples.append(sample)
            if sample_elements not in treated_condition:
                treated_condition.append(sample_elements)

    # Create matching matrix
    matrix = np.zeros((len(controls_condition), len(treated_condition)))

    # Checking the similarity between controls and treatments
    for i, control_sample in enumerate(controls_condition):
        for j, treated_sample in enumerate(treated_condition):
            for k, column in enumerate(columns):
                # If the level is not defined, continue or it is Replicate
                if '.0' in control_sample[k] or column == 'Replicate':
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

    # Grouping controls and treatments by subproject
    groups_list = [[] for x in range(len(controls_condition))]
    done = []
    for i in range(len(matrix)):
        max_v = max(matrix[i])
        # Add control samples
        for control in control_samples:
            if control in done:
                continue
            control_elements = control.split('_')
            del control_elements[2]
            if control_elements == controls_condition[i]:
                done.append(control)
                groups_list[i].append(control)
        
        # Add treated samples
        for j in range(len(matrix[i])):
            if matrix[i][j] == max_v:
                for treated in treated_samples:
                    if treated in done:
                        continue
                    treated_elements = treated.split('_')
                    del treated_elements[2]
                    if treated_elements == treated_condition[j]:
                        groups_list[i].append(treated)
                        done.append(treated)
    
    return groups_list


## MAIN PROGRAM

def main():
    '''
    Main program
    '''
    
    # Get and check arguments
    parser = argparse.ArgumentParser(prog='Divide_filter_by_subproject.py', 
                                     description='''This program divides the project count tables by
                                     subproject writing the table section of each subproject in a
                                     separate file. To do this, it uses the metadata of each project to
                                     group the columns of the treatment samples with the columns of the
                                     control samples with which they are associated. Finally, it extracts
                                     each group of columns from the original table (previously inserted into
                                     an SQLite database) and writes them to a new separate "CSV" file while
                                     selecting those sequences that have 5 or more counts in at least 5
                                     samples. In case the subproject has less than 5 samples, all the subproject
                                     samples must have at least 5 counts.''',  
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-c', '--path-counts', type=str, nargs=1,
                        help='Asbolute path of the directory where the project counts tables files are located')
    parser.add_argument('-m', '--metadata', type=str, nargs=1,
                        help='Asbolute path of the directory where the metadata tables files are located.')
    parser.add_argument('-r', '--path-results', type=str, nargs=1,
                        help='Asbolute path of the directory in which to save the results.')
    args = parser.parse_args()

    try:
        path_project_counts = args.path_counts[0]
        metadata = args.metadata[0]
        path_out = args.path_results[0]
    except:
        print('ERROR: You have inserted a wrong parameter or you are missing a parameter.')
        parser.print_help()
        sys.exit()
    

    # Get the project and the species name
    project = os.path.basename(path_project_counts).replace("_fusionCounts_RF", "")
    species = os.path.basename(os.path.dirname(path_project_counts)).replace("_fusionCounts_RF", "")

    # Create temporary directory
    os.system(f'mkdir -p ./tmp_{project}')

    # Get project name and path
    path_table = f'{path_project_counts}/fusion_abs-outer.csv'
                       
    # Table_name. Replace "-" with "_". SQLite does not accept the character "-".
    table_name = project.replace("-", "_")

    # Iterate paths list
    print(f'Separating {project} ({species}) project into subprojects...')

    # Get file columns
    with open(path_table, 'r') as file:
        columns = file.readline().rstrip().split(",")

    # Insert counts table into SQLite database
    id_dic = insert_to_database(f'./tmp_{project}/{project}_database.db', table_name, path_table, columns)

    # Obtain subprojects groups
    project_metadata = f'{metadata}/{species}_m_{project}.txt'
    groups_list = divide_project_by_subproject(columns, project_metadata)
    
    # Create project path out 
    project_path_out = f'{path_out}/{species}/{project}'
    os.system(f'mkdir -p {project_path_out}')

    # Write the subprojects in different files
    for i, group in enumerate(groups_list):
        file_path_out = f'{project_path_out}/{project}_{str(i + 1)}.csv'
        sql_to_csv(f'./tmp_{project}/{project}_database.db', table_name, file_path_out, group, id_dic)

    # Delete SQlite database
    os.system(f'rm -f -r ./tmp_{project}')
    print(f'{project} ({species}) done!')

## CALL THE MAIN PROGRAM

if __name__ == '__main__':
    main()
