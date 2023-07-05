#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#******************************************************************************
#  
#   08-Build_presence_absence_table.py
#
#   This program creates a table representing the presence (1) or absence
#   (0) of the miRNA families annotated in the different stress events.
#   The stress events are placed in the columns using identifiers created
#   by the program, where the species, the stress and the event in question
#   are considered. Meanwhile, the rows show the miRNA families annotated
#   in all the experiments considered. This table indicates that at least
#   one of the miRNA family members has been differentially expressed in a
#   given stress event with a 1,while the absence of these miRNAs is
#   indicated by a 0. This program generates two files, one of them being
#   the presence-absence table and the other a table where the identifiers
#   created from the program are related to the stress events to which they
#   are related, among other things.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 13/02/2023
#
#******************************************************************************


## IMPORT MODULES

import argparse
from natsort import natsorted
import os
import pandas as pd
import sys
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


## FUNCTIONS

def GetEventsIds(path_in: str, path_out: str):
    '''
    This function creates an identifier for each stress event. These
    identifiers consist of 3 numbers separated by ".". The first one refers to
    the species, the second one to the type of stress and the third one to the
    stress event. Finally, this function returns a dictionary whose key is the
    name of the file associated to that stress event (e.g. PRJNA277424_1_2)
    and as value a list containing the new generated identifier, the species,
    the stress and the stress event itself.

    Parameters
    ----------
    path_in : str
        Absolute path of the differential expression analysis summary file,
        which relates the file name to the stress event. 
    path_out : str
        Absolute path of the output file of the table with the new identifiers,
        files name, species, stresses and stress events.

    Returns
    -------
    dict
        Dictionary with the new identifiers, the files name, species, stresses
        and stress events.
    '''

    # Read DEA summary table
    comp_table = pd.read_csv(path_in, usecols=[0,1,2,3])
    comp_table = comp_table.reset_index()  # make sure indexes pair with number of rows

    # Get species-id dictionary
    species_list = comp_table['species'].unique().tolist()
    dic_species = { sp : i + 1 for i, sp in enumerate(species_list) }

    # Get stress-id dictionary
    comp_array = comp_table['Comparison']
    stress_list  = sorted(set([c.split("_")[2] for c in comp_array]))
    dic_stress = { sts : j + 1 for j, sts in enumerate(stress_list) }

    # Get id-event dictionary
    dic_events = {}
    dic_count = {}
    for index, row in comp_table.iterrows():

        # Get species and stress
        sp = row['species']
        stress = row['Comparison'].split("_")[2]
        
        # Build species-stress id
        id = f'{str(dic_species[sp])}.{str(dic_stress[stress])}'
        
        # Build event id
        if id not in dic_count:
            dic_count[id] = 1
            event_id = f'{id}.{str(dic_count[id])}'
        else:
            dic_count[id] += 1
            event_id = f'{id}.{str(dic_count[id])}'
        
        # Get experiment file name
        exp_filename = f'{row["Project"]}_{"_".join(str(row["Experiment"]).split("."))}'

        # Save results in dictionary
        dic_events[exp_filename] = (event_id, sp, stress, row['Comparison'])

    # Save ids in new file
    with open(path_out, 'a') as file_out:
        # Write header
        file_out.write('File,ID,Species,Stress,Comparison\n')
        # Write body
        for event_id in dic_events:
            file_out.write(event_id + ',' + dic_events[event_id][0]+ ',' +
                           dic_events[event_id][1]+ ',' +
                           dic_events[event_id][2] + ',' +
                           dic_events[event_id][3] + '\n')
    
    return dic_events


def GetmiRNAFamilies(path_in_annot: str, dic_id_event: dict, non_significant: bool=False):
    '''
    This function receives as parameters the absolute path of the differentially
    expressed miRNAs annotation directory and a dictionary in which the name of the
    file associated with a given stress event (e.g. PRJNA277424_1_2) is related
    to its new identifier (1.4.7), the species, the stress and the stress event
    itself. From these two variables, this function creates a dictionary with
    the miRNA families (values) that are differentially expressed in each
    stress event (key) and a list with all the annotated miRNA families.

    Parameters
    ----------
    path_in_annot : str
        Absolute path of the differentially expressed miRNAs annotation directory.
    dic_id_event : dict
        Dictionary in which the name of the file associated with a given stress
        event (key) is related to its new identifier (1.4.7), the species, the
        stress and the stress event (Values stored in a list in the same order).
    non_significant : bool
        Boolean specifying whether working with miRNAs that are not
        differentially expressed significantly. (optional. Default is False)

    Returns
    -------
    dict
        Dictionary with the miRNA families (values) that are differentially
        expressed in each stress event (key)
    list
        List with all the annotated miRNA families.
    '''

    # Iterate species
    miRNAs = []
    exp_miRNAs = {}
    species_list = os.listdir(path_in_annot)
    for species in species_list:

        # Discard summary file
        if species != 'summary.csv' and species != 'README.txt':

            # Iterate projects
            projects_list = os.listdir(f'{path_in_annot}/{species}')
            for project in projects_list:

                # Create project path
                project_path = f'{path_in_annot}/{species}/{project}/01-DEA_results_annot'

                # If working with non-significant miRNAs...
                if non_significant:

                    # Create project path for significantly differentially expressed miRNAs
                    path_sig_annot, data_dir = os.path.split(path_in_annot)
                    data_dir = data_dir.replace("nosig", "sig")
                    project_path_sig = f'{path_sig_annot}/{data_dir}/{species}/{project}/01-DEA_results_annot'

                # Iterate files
                files_list = os.listdir(project_path)
                for file in files_list:

                    # Get the experiment name
                    exp_name = file.rstrip(".csv")

                    # Get miRNA families of the experiment
                    miRNAs_temp = pd.read_csv(f'{project_path}/{file}', usecols=[13])['general_annot'].tolist()

                    # If working with non-significant miRNAs, remove from the list all miRNAs that are significant.
                    if non_significant:

                        # Check if the path exists
                        if os.path.exists(f'{project_path_sig}/{file}'):

                            # Get list of significant miRNA families
                            table_sig = pd.read_csv(f'{project_path_sig}/{file}', usecols=[13])
                            list_sig = list(set(table_sig.iloc[:, 0].values.tolist()))
                        
                            # Discard significant miRNA families
                            miRNAs_temp = [miRNA for miRNA in miRNAs_temp if miRNA not in list_sig]

                    # Get unique miRNAs
                    miRNAs += set(miRNAs_temp)

                    # Get the event id of the experiment
                    id = dic_id_event[exp_name][0]
                    exp_miRNAs[id] = set(miRNAs_temp)

                # Obtain unique miRNA families
                miRNAs = sorted(list(set(miRNAs)))

    return exp_miRNAs, miRNAs


def CreatePresenceAbsenceTable(dic_miRNAs_event: dict, miRNAs_list: list, path_out: str):
    '''
    This program creates a table representing the presence (1) or absence (0)
    of the miRNA families annotated in the different stress events. The stress
    events are placed in the columns using identifiers created by the program,
    where the species, the stress and the event in question are considered.
    Meanwhile, the rows show the miRNA families annotated in all the experiments
    considered. This table indicates that at least one of the miRNA family
    members has been differentially expressed in a given stress event with a 1,
    while the absence of these miRNAs is indicated by a 0.

    Parameters
    ----------
    dic_miRNAs_event : dict
        Dictionary with the miRNA families (values) that are differentially
        expressed in each stress event (key)
    miRNAs_list : list
        List with all the annotated miRNA families.
    path_out : str
        Absolute path of the file in which the presence-absence table will
        be stored.
    '''

    # Sort miRNAs and ids alphabetically
    ids_list = natsorted(list(dic_miRNAs_event.keys()))
    miRNAs_list = natsorted(miRNAs_list)

    # Iterate events ids
    final_table = pd.DataFrame()
    for id in ids_list:
        # Iterate miRNAs
        values_list = []
        for miRNA in miRNAs_list:
            if miRNA in dic_miRNAs_event[id]:
                values_list.append(1)
            else:
                values_list.append(0)

        # Add event column to final table
        final_table[id] = values_list
    
    # Add miRNA families as row names
    final_table.index = miRNAs_list
    final_table.index.name = 'miRNA_fam'

    # Add total row and column
    #final_table.loc['Total'] = final_table.sum(numeric_only=True, axis=0)
    #final_table.loc[:,'Total'] = final_table.sum(numeric_only=True, axis=1)

    # Add number of events and miRNA families
    #num_miRNAs = len(final_table) - 1 # Discard colnames (-1)
    #num_events = len(final_table.columns)
    #final_table.loc['Num_miRNAs'] = [num_miRNAs for i in range(num_events)]
    #final_table.loc[:,'Num_miRNAs'] = [num_events for i in range(num_miRNAs + 2)]
    
    # Save table in csv file
    final_table.to_csv(path_out, sep=",")


## MAIN PROGRAM

def main():
    '''
    Main program
    '''
    
    # Arguments
    parser = argparse.ArgumentParser(prog='Build presence-absence table', 
                                     description=''' This program creates a table representing the \
                                        presence (1) or absence (0) of the miRNA families annotated \
                                        in the different stress events. The stress events are placed \
                                        in the columns using identifiers created by the program, where \
                                        the species, the stress and the event in question are considered. \
                                        Meanwhile, the rows show the miRNA families annotated in all the \
                                        experiments considered. This table indicates that at least one of \
                                        the miRNA family members has been differentially expressed in a \
                                        given stress event with a 1, while the absence of these miRNAs is \
                                        indicated by a 0. This program generates two files, one of them \
                                        being the presence-absence table and the other a table where the \
                                        identifiers created from the program are related to the stress \
                                        events to which they are related, among other things.''',  
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-d', '--path-dea', type=str, nargs=1,  
                        help='Absolute path of the differential expression \
                        analysis summary file, which relates the file name \
                        to the stress event. ')
    parser.add_argument('-a', '--path-annot', type=str, nargs=1, 
                        help='Absolute path of the differentially expressed \
                        miRNAs annotation directory.')
    parser.add_argument('-o', '--path-dir-out', type=str, nargs=1,
                        help='Absolute path of results directory.')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    # Check the arguments and store them in a variable
    try:
        path_in_sum_dea = args.path_dea[0]
        path_in_annot = args.path_annot[0]
        path_out = args.path_dir_out[0]

    except:
        print('ERROR: You have inserted a wrong parameter or you are missing a parameter.')
        parser.print_help()
        sys.exit()

    # Create output directory
    os.makedirs(path_out)

    # Get ids-event dictionary
    print('Creanting stress event ids...')
    dic_id_event = GetEventsIds(path_in_sum_dea, f'{path_out}/ids_table.csv')
    print('Done!')

    # Iterate data list (sig or nosig)
    data_list = os.listdir(path_in_annot)
    for data in data_list:
        
        # Get data path and suffix (sig or nosig)
        path_annot_data = f'{path_in_annot}/{data}'
        suffix = data.split("_")[-1]

        # Significant or non-significant data?
        if suffix == "nosig":
            non_sig =  True
        else:
            non_sig =  False
        
        # Get miRNAs and the events in which they are differentially expressed
        print('Relating stress event identifiers to miRNA families...')
        dic_miRNAs_event, miRNAs_list = GetmiRNAFamilies(path_annot_data, dic_id_event, non_sig)
        print('Done!')

        # Create Presence-Absence table
        print('Creating presence-absence table...')
        CreatePresenceAbsenceTable(dic_miRNAs_event, miRNAs_list, f'{path_out}/presence_absence_table_{suffix}.csv')
        print('Done!')


## CALL THE MAIN PROGRAM

if __name__ == '__main__':
    '''
    Call to the main program
    '''
    main()
