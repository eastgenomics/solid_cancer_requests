import concurrent.futures
import dxpy as dx
import pandas as pd
import sys


def get_002_TSO500_projects_in_period(first_date, last_date):
    """
    Gets all the 002 TSO500 projects between the period given
    Parameters
    ----------
    first_date : str
        e.g. '2023-07-01'
    last_date : str
        e.g. '2023-11-01'
    Returns
    -------
    tso_projs_ids : list
        list of project IDs within the period
    """
    #Search projs in the period starting with 002 and ending with TSO500
    #Return only the ID field from describe
    tso_projs_response = list(dx.find_projects(
        level='VIEW',
        created_before=last_date,
        created_after=first_date,
        name="002*TSO500",
        name_mode="glob"
    ))

    #Get only the IDs from each response dict
    tso_projs_ids = [proj['id'] for proj in tso_projs_response]

    return tso_projs_ids


def get_name_of_TSO_folder(project_id):
    """
    Get the name of the main analysis folder in the TSO500 project
    to enable searching a specific folder for the CombinedVariantOutput
    files later

    Parameters
    ----------
    project_id : str
        ID of the project being searched

    Returns
    -------
    folder_name : str
        name of the folder within the output folder e.g. 'TSO500-230724_0472"
    """
    folder_names = list(
        dx.bindings.dxfile_functions.list_subfolders(
            project=project_id,
            path='/output/',
            recurse=False
        )
    )
    #Check if there is only one run in the project and get that folder name
    folder_name = []
    #in case there is only one run get that folder name (without path)
    if len(folder_names) == 1:
        folder_name = folder_names[0].removeprefix('/output/')
    #if more get the one with reports workflow
    if len(folder_names) > 1:
        print(f"Warning: more than one folder found for project {project_id}")
        for folder in folder_names:
            #list subfolders in each project
            sub_folders = list(
                dx.bindings.dxfile_functions.list_subfolders(
                project=project_id,
                path=folder,
                recurse=False
                )
            ) 
            #check for all the runs and get the one with the reports folder
            for sub_folder in sub_folders:
                if "TSO500_reports_" in sub_folder:
                    #get only the folder with date name
                    folder_name = sub_folder.split('/')[2]

    return folder_name


def get_combinedvariantoutput_files(project_id, folder_name):
    """
    Find the CombinedVariantOutput files in the relevant path for each project,
    subsetting to only files with 8471 and 8475 in the name as they are DNA

    Parameters
    ----------
    project_id : str
        the proj ID as a string
    Returns
    -------
    combined_var_files : list
        list of dicts containing info about the files found
    """
    #List all of the CombinedVariantOutput files in the specific folder
    #of that project. Only find files with 8471 or 8475 in the name (DNA)
    combinedvar_files_response = list(
        dx.find_data_objects(
            project=project_id,
            folder=(
                f"/output/{folder_name}/eggd_tso500/"
            ),
            recurse=True,
            name=".*[8471|8475]_CombinedVariantOutput.tsv",
            name_mode='regexp',
            classname='file',
            describe={
                'fields': {
                    'name': True,
                    'archivalState': True
                }
            }
        )
    )
    #Get only these fields in a simple dict for easier querying later
    #and to make the merged list of the hundreds of files slightly
    #smaller later
    combined_var_files = [
        {
            'project': x['project'],
            'id': x['id'],
            'name': x['describe']['name'],
            'archive_state': x['describe']['archivalState']
        } for x in combinedvar_files_response
    ]

    return combined_var_files


def find_archived_files(combined_var_files):
    """
    Check the archivalState of all of the CombinedVariantOutput files
    and return any files which are not live

    Parameters
    ----------
    combined_var_files : list
        list of dicts containing info about each file

    Returns
    -------
    archived_files : list
        list of file IDs of files which are not live
    """
    archived_files = []
    for x in combined_var_files:
        if x['archive_state'] != 'live':
            archived_files.append({'project': x['project'], 'id': x['id']})

    return archived_files


def unarchive_files(archived_files):
    """
    Call unarchive on any files which are not live

    Parameters
    ----------
    archived_files : list
        list of file IDs of non-live files
    """
    for idx, file in enumerate(archived_files):
        print(f"Checking file {file}{idx + 1}/{len(archived_files)}")
        file_object = dx.DXFile(file.get('id'), project=file.get('project'))
        file_object.unarchive()


def read_to_df(file_dict):
    """
    Read from DNA Sample ID to Gene Amplifications 
    (excluding) part of the CombinedVariantOutput tsv
    (from a DNAnexus file ID) into a pandas dataframe

    Parameters
    ----------
    file_dict : dict
        dict containing information about the CombinedVariantOutput file

    Returns
    -------
    msi_tmb : pd.DataFrame
        pandas df of metrics and value for the sample
    """
    file_id = file_dict['id']
    
    with dx.open_dxfile(file_id, mode='r') as dx_file:
        # Read TSV to pandas df
        data = pd.read_csv(
            dx_file, sep='\t', header=None, names=list(range(11))
        )
        # Get indexes of section of TSV we want
        start_string = 'DNA Sample ID'
        idx_start = data.index[data[0] == start_string][0]
        end_string = '[Gene Amplifications]'
        idx_end = data.index[data[0] == end_string][0]

        # Subset df to just those columns and first two columns
        msi_tmb = data.loc[idx_start : idx_end-1][[0, 1]]
        # Change column names to the Metric and value associated
        msi_tmb.columns = ['Metric', 'Value']
        #insert random column numbered
        #to make transpose easier as columns need different names
        #remove later
        msi_tmb.insert(0, 'to_remove', range(1, 1 + len(msi_tmb)))
        #create new df with transposed data
        new = msi_tmb.set_index('to_remove').T

    return new


def concurrent_read_tsv(list_of_file_dicts, workers):
    """
    Concurrently read in the info selected to a df for each
    TSV file

    Parameters
    ----------
    list_of_file_dicts : list
        list of dicts, each dict containing info on a 
        CombinedVariantOutput file
    workers : int
        Number of workers

    Returns
    -------
    single_df : pd.DataFrame
        single df with the gene 
        amplification data from all samples
    """
    list_of_dfs = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        concurrent_jobs = {
            executor.submit(
                read_to_df, file_dict
            ): file_dict for file_dict in list_of_file_dicts
        }
        for future in concurrent.futures.as_completed(concurrent_jobs):
            try:
                data = future.result()
                list_of_dfs.append(data)
            except Exception as exc:
                print(
                    f"Error reading file for {concurrent_jobs[future]}: {exc}"
                )
    #Concat all the dfs in the list to one
    single_df = pd.concat(list_of_dfs)

    return single_df


def main():
    #get tsoo500 projects
    tso_projs_ids = get_002_TSO500_projects_in_period("2023-07-15","2024-05-15")

    #list of file dict
    combined_var_files = []
    for proj in tso_projs_ids:
        folder_to_search = get_name_of_TSO_folder(proj)
        combined_var_files.extend(get_combinedvariantoutput_files(
        proj, folder_to_search
    ))
        
    #find files that are archived
    archived_files = find_archived_files(combined_var_files)

    #print archived files, unarchive them and get status
    if archived_files:
        print(
            f"{len(archived_files)}files are still archived. Calling "
            "unarchive on them and exiting. Please re-run later when files are"
            "unarchived"
        )
        unarchive_files(archived_files)
        sys.exit(1)
    else:
        print("All files are unarchived. Reading into dataframes and concatenating")

    #create df with selected info from all files
    all_samples_gene_df = concurrent_read_tsv(combined_var_files, 12)
    #drop duplicates from transposing - keep first to act as column
    all_samples_gene_df = all_samples_gene_df.drop_duplicates(keep='first')
    #get rid of metrics that are not important
    all_samples_gene_df = all_samples_gene_df.iloc[:, [0,8,20,21,22,25,26,27]] 
    #drop rows that are empty
    all_samples_gene_df.dropna(inplace = True) 
    #save df into csv file
    all_samples_gene_df.to_csv('TSO500_TMB_MSI_all.csv', index=False, header=False)

        #Get final matches between helpdesk samples and all obtained
    #get samples from helpdesk and name column
    #copied from helpdesk to a excel - removed SP with sed 
    helpdesk_table = pd.read_csv('helpdesklist_nosp.csv', sep=',')
    helpdesk_table.columns = (['sample'])

    #read samples
    all_samples = pd.read_csv('TSO500_TMB_MSI_all.csv', sep=',')

    #convert sample name
    all_samples['sample_name'] = all_samples["DNA Sample ID"].str.split('-').str[1]

    #Initialize an empty list to store matching sample names
    matching_samples = []

    #Iterate over each sample_name in all_samples
    for i in all_samples['sample_name']:
        #Iterate over each sample in helpdesk_table
        for l in helpdesk_table['sample']:
            #Convert l and i to strings before stripping
            l_str = str(l)
            i_str = str(i)
            #Check if the stripped versions of l and i are equal
            if l_str.strip() == i_str.strip():
                #If a match is found, append the sample name 
                # to the list and break the loop
                matching_samples.append(i)
                print('true')
                break  
        else:
        
            print('')  

    #Filter all_samples based on matching sample names
    results = all_samples[all_samples["sample_name"].isin(matching_samples)]
    #save matches and all to csv
    results.iloc[:,1:].to_csv('TSO500_TMB_MSI_matches.csv', index=False)
    all_samples.iloc[:,1:].to_csv('TSO500_TMB_MSI_final.csv', index=False)
        
if __name__ == "__main__":
    main()


