#!/usr/bin/env python
# coding: utf-8

# In[1]:


import mapping_tools as mt
import numpy as np
import sys
import os
import gzip
import re
import pandas as pd
from timeit import default_timer as timer
from VEPcrossref import VEPfileCrossrefGenerator as cr


# In[ ]:


# INTERFACES PART


# In[2]:


# Extract the info corresponding to the prot ID (Interface parse)

#sys.argv[1]
#sys.argv[2]

protID = 'ENSP00000482258'
int_db_file = '/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/raw_data/interfaces_mapped_to_v94.csv'

  
interfacesDB = open(int_db_file, "r")
cols = interfacesDB.readline().strip().split(" ")
#prot_interface = mt.ProtIntDB_parser(interfacesDB, protID)
prot_interface = mt.parser(interfacesDB,  protID, cols, " " )


# In[10]:


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


# In[ ]:


####################### VARIANTS PART #####################################################


# In[32]:


# Translate protID to gene ID
#subset_prot_interface


# In[4]:


biomartdb = open("/home/vruizser/PhD/2018-2019/git/PDBmapper/project/gene_transcript_protein_ens_ids.txt", "r")
ensemblIDs = mt.ensemblID_translator(biomartdb, protID)


# In[5]:


# get the corresponding VEP file
crossref_file = open("/home/vruizser/PhD/2018-2019/git/PDBmapper/project/geneids_VEPfiles_crossref.txt", "r")
VEP_filename = mt.VEP_getter(crossref_file, ensemblIDs["geneID"])


# In[6]:


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

 

    
    
    return


##########################
# execute main function  #
##########################
if __name__ == '__main__':
    main()

