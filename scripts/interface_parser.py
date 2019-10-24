# Import necesary modules
import numpy as np
import sys
import os
import re
import pandas as pd
from scripts.explode import explode


# Extract the info corresponding to the prot ID (Interface parse)
def reshape(annoint):
    '''Parse input interfaces database to put it in the right format.
    Parameters
    ----------
    protID : str
        Ensemble protein id 
    interfacesDB_filepath : str
        DESCRIPTION MISSING!!
    Returns
    -------
    subset_interfaces_db
        DESCRIPTION MISSING!!
    '''
    # store subspace
    sub_annoint = annoint[['pdb.id',
                           'ensembl.prot.id',
                           'temp.chain',
                           'int.chain',
                           'interaction',
                           'resid_sseq',
                           'mapped.real.pos',
                           'pdb.pos']]
   
    # put it into right format
    sub_annoint.columns = \
        sub_annoint.columns.str.replace('\\.', '_')

    for col in ('resid_sseq', 'mapped_real_pos', 'pdb_pos'):
        sub_annoint.loc[:, col] = \
            sub_annoint[col].str.replace('-', ',')
        sub_annoint.loc[:, col] = \
            sub_annoint[col].str.split(',')

    sub_annoint = explode(sub_annoint,
                          ['resid_sseq',
                           'mapped_real_pos',
                            'pdb_pos'])
    sub_annoint.rename(columns={'mapped_real_pos':
                                          'Protein_position'}, inplace=True)
    # create region id for setID file
    sub_annoint['region_id'] = sub_annoint['pdb_id'] + '_' + \
        sub_annoint['ensembl_prot_id'] + '_' + \
        sub_annoint['temp_chain'] + '_' + \
        sub_annoint['int_chain'] + '_' + \
        sub_annoint['interaction']

    return sub_annoint
