#!/usr/bin/env python
# coding: utf-8

# Import necesary modules 
import mapping_tools as mt
import numpy as np
import sys
import os
import gzip
import re
import pandas as pd
from timeit import default_timer as timer
from VEPcrossref import VEPfileCrossrefGenerator as cr

# Extract the info corresponding to the prot ID (Interface parse)

#protID = 'ENSP00000482258'
#int_db_file = '/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/raw_data/interfaces_mapped_to_v94.csv'


def PDBmapper(protID, interfacesDB_filepath):
    """Generate setID file.
    
    Parameters
    ----------
    protID : str
        input path of the VEP file
    interfacesDB_filepath : str
        chosen name of the output file

    Returns
    -------
    setID.File
        write data frame to a txt file with two columns: one is the gene ids and the other one the VEP file
    data.frame 
        kasklkssklasdkmaskldmakmdalkmd
    """

    # read interfaces file
    interfacesDB = open(interfacesDB_filepath, "r")
    # get colnames of interfaces file
    intDB_colnames = interfacesDB.readline().strip().split(" ")
    # parse information related to protein ID
    prot_interface = mt.parser(interfacesDB,  protID, intDB_colnames, " " )
    # store 
    subset_prot_interface = prot_interface[['pdb.id',
                                            'ensembl.prot.id',
                                            'temp.chain',
                                            'int.chain',
                                            'interaction',
                                            'resid_sseq',
                                            'mapped.real.pos',
                                            'pdb.pos' ]]

    # put it into right format
    subset_prot_interface.columns = subset_prot_interface.columns.str.replace("\\.", "_")

    for col in ("resid_sseq", "mapped_real_pos", "pdb_pos"):
        subset_prot_interface.loc[ : , col] = subset_prot_interface[col].str.replace("-", ',')
        subset_prot_interface.loc[ : , col] = subset_prot_interface[col].str.split(',')

    subset_prot_interface = mt.explode(subset_prot_interface, ["resid_sseq", "mapped_real_pos", "pdb_pos"])

    subset_prot_interface.rename(columns={'mapped_real_pos':'Protein_position'}, inplace=True)

    subset_prot_interface['region_id'] = subset_prot_interface['pdb_id'] + '_' + subset_prot_interface['ensembl_prot_id'] + '_' + subset_prot_interface['temp_chain'] + '_' + subset_prot_interface['int_chain'] + '_' + subset_prot_interface['interaction']

    biomartdb = open("/home/vruizser/PhD/2018-2019/git/PDBmapper/project/gene_transcript_protein_ens_ids.txt", "r")
    ensemblIDs = mt.ensemblID_translator(biomartdb, protID)



    # get the corresponding VEP file
    crossref_file = open("/home/vruizser/PhD/2018-2019/git/PDBmapper/project/geneids_VEPfiles_crossref.txt", "r")
    VEP_filename = mt.VEP_getter(crossref_file, ensemblIDs["geneID"])



VEP_dir= "/home/vruizser/PhD/2018-2019/PanCancer_data/vep_output/"
VEPfile = open(VEP_dir + VEP_filename, "r")
geneID = ensemblIDs['geneID']
cols = pd.read_csv(VEPfile, nrows = 0, skiprows = 42, sep = "\t").columns
VEP = mt.parser(VEPfile,  geneID, cols, "\t" )


# In[9]:


setID = MapVariantToPDB(VEP, subset_prot_interface, 'region_id')


# In[ ]:


# define main function to execute the previous defined functions together
def main():

    protID = sys.argv[1]
    interfacesDB_filepath = sys.argv[2]

    PDBmapper(protID, interfacesDB_filepath)
    return


##########################
# execute main function  #
##########################
if __name__ == '__main__':
    main()

