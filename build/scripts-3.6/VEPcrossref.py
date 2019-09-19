#!/usr/bin/python3
import pandas as pd
import os

def VEPfileCrossrefGenerator(VEP_file_path, crossref_filename):
    """Generate file of genes IDs and the corresponding VEP file where you can find them.
    
    Parameters
    ----------
    VEP_file_path : str
        input path of the VEP file
    crossref_filename : str
        chosen name of the output file

    Returns
    -------
    data frame
        write data frame to a txt file with two columns: one is the gene ids and the other one the VEP file
    """

    #read the VEP file skipping 42 rows of information
    VEPfile = pd.read_csv(VEP_file_path, skiprows=42, error_bad_lines=False, sep='\t')
    #extract the unique gene ids 
    gene_unique_ids = pd.DataFrame(VEPfile.Gene.unique(), columns=['ids'])
    #create 
    gene_unique_ids['VEPfile'] = os.path.basename(VEP_file_path)
    #append results to new file 
    if os.path.isfile('.'+ crossref_filename) :

        gene_unique_ids.to_csv(crossref_filename, mode='a', sep='\t', header=False, index = False)

    else :

        gene_unique_ids.to_csv(crossref_filename, mode = 'a', sep='\t', header=True, index = False)
        
