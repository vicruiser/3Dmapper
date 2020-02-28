# coding: utf-8
# Import necesary modules
import numpy as np
import sys
import os
import re
import pandas as pd
from .explode import explode

# Classify location of input mutations


# Extract the info corresponding to the prot ID (Interface parse)
# def reshape(annoint):
# '''
# Parse input interfaces database to put it in the right format.

# Parameters
# ----------
# annoint : df
#     input interfaces file

# Returns
# -------
# sub_annoint
#     reshaped interfaces file
# '''
# # store subspace
# sub_annoint = annoint[['PDB_code',
#                        'Protein_accession',
#                        'PDB_chain',
#                        'PDB_interacting_chain',
#                        'Interaction',
#                        'Protein_aa',
#                        'Protein_pos',
#                        'PDB_pos',
#                        'Pident']]

# # put it into right format
# # sub_annoint.columns = \
# #     sub_annoint.columns.str.replace('\\.', '_')
# # if database only contains one row, explode is not needed
# if any(sub_annoint[['Protein_aa',
#                     'Protein_pos',
#                     'PDB_pos']].stack().str.contains("-", na=False)) is True:

#     for col in ['Protein_aa', 'Protein_pos', 'PDB_pos']:

#         sub_annoint.loc[:, col] = sub_annoint[col].str.replace(
#             '-', ',')
#         sub_annoint.loc[:, col] = sub_annoint[col].str.split(',')

#     sub_annoint = explode(
#         sub_annoint, ['Protein_aa', 'Protein_pos', 'PDB_pos'])
# elif 1 > 0:
#     sub_annoint = explode(
#         sub_annoint, ['Protein_aa', 'Protein_pos', 'PDB_pos'])
# else:
#     sub_annoint[['Protein_aa', 'Protein_pos', 'PDB_pos']] = \
#         sub_annoint[['Protein_aa', 'Protein_pos', 'PDB_pos']].astype(str)

# # # change name of columns
# # sub_annoint.rename(columns={'mapped_real_pos':
# #                             'Protein_position'}, inplace=True)

# # create region id for setID file
# # sub_annoint['interface_id'] = sub_annoint['pdb_id'].astype(str) + '_' + \
# #     sub_annoint['ensembl_prot_id'].astype(str) + '_' + \
# #     sub_annoint['temp_chain'].astype(str) + '_' + \
# #     sub_annoint['int_chain'].astype(str) + '_' + \
# #     sub_annoint['interaction'].astype(str)
# # drop unnecesary columns
# columns = ['pdb_id',
#            'ensembl_prot_id',
#            'temp_chain',
#            'int_chain',
#            'interaction']
# sub_annoint.drop(columns, inplace=True, axis=1)
# return sub_annoint
