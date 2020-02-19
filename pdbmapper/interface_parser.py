# coding: utf-8
# Import necesary modules
import numpy as np
import sys
import os
import re
import pandas as pd
from .explode import explode


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
    sub_annoint = annoint[['PDB_code',
                           'ensembl.prot.id',
                           'temp.chain',
                           'int.chain',
                           'interaction',
                           'resid_sseq',
                           'mapped.real.pos',
                           'pdb.pos',
                           'pident']]

    # put it into right format
    sub_annoint.columns = \
        sub_annoint.columns.str.replace('\\.', '_')
    # if database only contains one row, explode is not needed
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

    # change name of columns
    sub_annoint.rename(columns={'mapped_real_pos':
                                'Protein_position'}, inplace=True)

    # create region id for setID file
    sub_annoint['interface_id'] = sub_annoint['pdb_id'].astype(str) + '_' + \
        sub_annoint['ensembl_prot_id'].astype(str) + '_' + \
        sub_annoint['temp_chain'].astype(str) + '_' + \
        sub_annoint['int_chain'].astype(str) + '_' + \
        sub_annoint['interaction'].astype(str)
    # drop unnecesary columns
    columns = ['pdb_id',
               'ensembl_prot_id',
               'temp_chain',
               'int_chain',
               'interaction']
    sub_annoint.drop(columns, inplace=True, axis=1)
    return sub_annoint
