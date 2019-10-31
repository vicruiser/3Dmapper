# Import necesary modules
import numpy as np
import sys
import os
import re
import pandas as pd
from scripts.explode import explode


# Extract the info corresponding to the prot ID (Interface parse)
def reshape(annoint):
    '''
    Parse input interfaces database to put it in the right format.

    Parameters
    ----------
    annoint : df
        input interfaces file

    Returns
    -------
    sub_annoint
        reshaped interfaces file
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
    # sub_annoint[['mapped_real_pos', 'pdb_pos']] = sub_annoint[[
    #     'mapped_real_pos', 'pdb_pos']].astype(str)

    # put it into right format
    sub_annoint.columns = \
        sub_annoint.columns.str.replace('\\.', '_')

    # if sub_annoint[sub_annoint.stack().str.contains("-")] != '':
    if any(sub_annoint[['resid_sseq', 'mapped_real_pos', 'pdb_pos']].stack().str.contains("-", na=False)) is True:

        for col in ['resid_sseq', 'mapped_real_pos', 'pdb_pos']:

            sub_annoint.loc[:, col] = sub_annoint[col].str.replace(
                '-', ',')
            sub_annoint.loc[:, col] = sub_annoint[col].str.split(',')

        sub_annoint = explode(sub_annoint,
                              ['resid_sseq',
                               'mapped_real_pos',
                               'pdb_pos'])
    else:
        sub_annoint[['resid_sseq', 'mapped_real_pos', 'pdb_pos']] = sub_annoint[[
            'resid_sseq', 'mapped_real_pos', 'pdb_pos']].astype(str)

    sub_annoint.rename(columns={'mapped_real_pos':
                                'Protein_position'}, inplace=True)

    # create region id for setID file
    sub_annoint['region_id'] = sub_annoint['pdb_id'] + '_' + \
        sub_annoint['ensembl_prot_id'] + '_' + \
        sub_annoint['temp_chain'].astype(str) + '_' + \
        sub_annoint['int_chain'].astype(str) + '_' + \
        sub_annoint['interaction']

    return sub_annoint
