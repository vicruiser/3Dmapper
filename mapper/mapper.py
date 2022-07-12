# coding: utf-8

# Import necesary modules
import sys
import os
import re
import glob
import pandas as pd
import numpy as np
#import dask.dataframe as dd

from .db_parser import parser
from .decorator import tags
from .explode import explode
from .explode2 import explode2
from .logger import get_logger
from .writefile import writefile

def mapper(prot_id,  gene_id, transcript_id, psdb, vardb, out_dir, pident, evalue, isoform, APPRIS, consequence, loc, var_id=None, csv=False, hdf=False):
    print(prot_id)
    '''
    Map interfaces and genomic anntoated positions and returns a
    setID.File, necessary input for SKAT. Additionaly, it creates
    another file with detailed information regarding the maped areas.

    Parameters
    ----------
    prot_id : str
        Ensembl protein ID
    gene_id : str
        Translated Ensembl protein ID
    psdb : str
        Directory where to find interface database
    vardb : str
        Directory where to find positions database
    out_dir : str
        Output directory
    pident : int
        Thershold of sequence identity (percertage).

    Returns
    -------
    setID.File
        txt file containing a data frame two columns corresponding to the
        analyzed interface id and the corresponding annotated genomic positions.
    InterfacePositions
        Same as setID.File but with additional information describing the
        interfaces and the positions.
    '''
    # log file
    logger = get_logger(' 3dmapper', out_dir)
    # amino acids letters
    aa = ['I', 'M','T','N', 'K', 'S', 'R', 'L', 
     'P', 'H', 'Q', 'V','A', 'D',
     'E','G','F', 'Y', 'C', 'W','X']
    # parse positions corresponding to the selected protein ID
    try:
        annovars = parser(transcript_id, vardb)     
        if consequence is not None:
            annovars = annovars[annovars['Consequence'].astype(
                str).str.contains('|'.join(consequence))]
            logger.info('Filter of features = ' + str(consequence))

            # if filter returns an empty df, raise error
            if annovars.empty is True:
                logger.error(
                    'positions could not be filtered by feature type = ' + consequence)
                raise IOError()
        else:
            consequence = ['all']

        if isoform is None:
            isoform = ['all']

        # filter by position type if one or more selected
        if var_id is not None:
            if 'Existing_variation' in annovars.columns:
                annovars = annovars[
                    (annovars.Uploaded_variation.astype(str) == str(var_id)) |
                    (annovars.Existing_variation.astype(str) == str(var_id))]
            else:
                annovars = annovars[
                    (annovars.Uploaded_variation.astype(str) == str(var_id))]
            logger.info('position \'' + str(var_id) + '\' has been selected.')
            # if filter returns an empty df, raise error
            if annovars.empty:
                logger.error(
                    'positions could not be filtered by position id \'' + str(var_id) + '\'')
                raise IOError()
        # for positions with high impact affecting several aminoacidic positions,
        # the protein position is a range. split the range to have each position
        # individually
        if any(annovars['Protein_position'].astype(str).str.contains(r'[0-9]-[0-9]')):
            # subset hight impact positions
            sub_df = annovars[annovars['Protein_position'].astype(str).str.contains(
                r'[0-9]-[0-9]')]
            # subset the remaining positions to concatenate afterwards
            remaining_df = annovars.drop(sub_df.index)
            # split the range or interval
            sub_df[['start', 'end']] = sub_df['Protein_position'].str.split(
                '-', expand=True)
            # sometimes the start or the end position of the interval is a
            # question mark. In that case, we take into account the
            # remaining value of the interval
            if any(sub_df['start'].str.contains('\?')):
                sub_df['start'] = np.where(sub_df['start'] == '?', sub_df['end'],
                                           sub_df['start'])
            if any(sub_df['end'].str.contains('\?')):
                sub_df['end'] = np.where(sub_df['end'] == '?', sub_df['start'],
                                         sub_df['end'])
            # create the range of numbers defined by the interval
            sub_df['Protein_position'] = sub_df.apply(lambda x: list(
                range(int(x['start']), int(x['end'])+1)), 1)
            # spread each individual position into one row
            sub_df = explode2(sub_df, ['Protein_position'])
            # drop unnecesary columns
            sub_df.drop(['start', 'end'], inplace=True, axis=1)
            # concatenate final result
            annovars = pd.concat([remaining_df, sub_df], sort=False)
            annovars = annovars.reset_index(drop=True)
    except IOError:
        annovars = False
     # parse interfaces corresponding to the selected protein ID
    try:
        psdf = parser(prot_id, psdb)
        if 'Pident' not in list(psdf.columns):
            logger.error(' Wrong structural data format. Header is missing')
            raise IOError()
        psdf = psdf.loc[psdf.Pident != 'Pident']
        logger.info('Protein features file of ' + prot_id + ' parsed.')
        # Get col PDB_position
        # if database is compacted and explode is needed
        columns_type = psdf.applymap(lambda x: isinstance(x, list)).all()
        columns_list = columns_type.index[columns_type].tolist()

        # cols_stack
        #cols_stack = psdf.apply(lambda x: x.astype(
        #    str).str.match(r'[a-zA-Z0.-9]+/[a-zA-Z0.-9]+'))
        #colsnames_stack = psdf.columns[cols_stack.any()].tolist()
        # add column for chimera script
        #if 'PDB_interacting_3D_position' in colsnames_stack:
            #psdf['Interface_interacting_positions'] = psdf['PDB_interacting_position']  
        #    psdf['Chimera_interacting_position'] = psdf['PDB_interacting_3D_position']
        #if 'Evalue' in colsnames_stack:
        #    colsnames_stack.remove('Evalue')
        #if 'PDB_code' in colsnames_stack:
        #    colsnames_stack.remove('PDB_code')
        #if 'PDB_3D_position' in colsnames_stack:
            #colsnames_stack.remove('PDB_3D_position')
        #    psdf['Chimera_3D_position'] = psdf['PDB_3D_position']
        #if 'Structure_feature_id' in colsnames_stack:
        #    colsnames_stack.remove('Structure_feature_id')
        #if any(colsnames_stack):
        #    psdf = explode(psdf , colsnames_stack, '/')
        #if any(columns_list):
        #    psdf = explode(
        #        psdf, columns_list, '/')
        #else:
        #    psdf[colsnames_stack] = \
        #        psdf[colsnames_stack].astype(str)
        # if default database is used minor modifications are needed
        if pident is not None:
            logger.info('Filtering interfaces by pident = ' +
                        str(pident) + '%.')
            # filter by pident
            pident = float(pident)  # from str to int
            psdf = psdf.loc[psdf.Pident.astype(float) >= pident]
            # if pident threshold is to high, the next maximum value of pident is
            # notified in log file
            if psdf.empty:
                alt_pident = psdf.loc[:, "Pident"].max()
                logger.error('Warning: for prot_id ' + str(pident) +
                             ', the variable "Pident" equal to ' +
                             str(pident) + ' is too high.\n A threshold lower than or equal to ' +
                             str(alt_pident) + ' would retrieve results.')

                raise IOError()
        if evalue is not None:
            logger.info('Filtering interfaces by evalue = ' +
                        str(evalue) + '%.')
            # filter by pident
            evalue = float(evalue)  # from str to int
            psdf = psdf.loc[psdf.Evalue.astype(float) >= evalue]
            # if pident threshold is to high, the next maximum value of pident is
            # notified in log file
            if psdf.empty:
                alt_evalue = psdf.loc[:, "Evalue"].min()
                logger.error('Warning: for prot_id ' + str(evalue) +
                             ', the variable "Evalue" equal to ' +
                             str(evalue) + ' is too low.\n A threshold higher than or equal to ' +
                             str(alt_evalue) + ' would retrieve results.')

                raise IOError()       
    except IOError:
        psdf = False

    if psdf is False and annovars is False:
        logger.error('Protein ' +
                     prot_id + 'could not be parsed.')
        raise IOError

    elif psdf is not False and annovars is False:
        logger.error('Protein ' +
                     prot_id + 'has no mapping positions.')
        raise IOError

    elif psdf is False and annovars is not False:

        if loc:
            # remove non protein coding positions
            annovars.drop_duplicates(inplace=True)
            try:
                coding_positions_index = annovars.Amino_acids.str.contains(
                    '|'.join(aa), regex=True, na=False)
                noncoding_positions = annovars.loc[~coding_positions_index]
            except:
                noncoding_positions = False
            # non-protein coding mutations
            if noncoding_positions is not False:
                if APPRIS is not None: 
                    noncoding_positions['APPRIS_isoform'] = APPRIS
                noncoding_positions['Mapping_position'] = 'Noncoding'
                writefile(prot_id, out_dir, float(pident), isoform, consequence,
                          noncoding_positions, 'NoncodingPositions', csv, hdf)
                unmapped_positions = annovars.loc[coding_positions_index]
            else:
                unmapped_positions = annovars

            if unmapped_positions.empty is False:
                unmapped_positions.drop_duplicates(inplace=True)
                if APPRIS is not None: 
                    unmapped_positions['APPRIS_isoform'] = APPRIS
                
                unmapped_positions['Mapping_position'] = 'Unmapped'
                writefile(prot_id, out_dir, float(pident), isoform, consequence,
                          unmapped_positions, 'UnmappedPositions', csv, hdf)
   
    elif psdf is not False and annovars is not False:
        # for sucessful merge, Protein_position column must be str type
        psdf['Protein_position'] = psdf['Protein_position'].astype(str)
        annovars['Protein_position'] = annovars['Protein_position'].astype(str)
        if APPRIS is not None:
            annovars['APPRIS_isoform'] = APPRIS
       # annovars['Uniprot_accession'] = Uniprot_id
        # Merge them both files
        mapped_positions = annovars.join(
            psdf.set_index('Protein_position'), on='Protein_position', how='inner')
        if isoform is None:
            isoform = ['all']
        if consequence is None:
            consequence = ['all']
        ###########################################################################
        # Locate rest of positions (mapping to a structure or not)
        ###########################################################################
        if loc:
            # remove already mapped positions
            try:
                left_positions = annovars.drop(list(set(mapped_positions.index)))
            except:
                left_positions = annovars

            # remove non protein coding positions
            if left_positions.empty is False:
                left_positions.drop_duplicates(inplace=True)
                coding_positions_index = left_positions.Amino_acids.str.contains(
                    '|'.join(aa), regex=True, na=False)
                noncoding_positions = left_positions.loc[~coding_positions_index]
                # non-protein coding mutations
                if noncoding_positions.empty is False:
                   # noncoding_positions = noncoding_positions.drop(
                   #     columns=['Uniprot_accession'])  # not needed
                    noncoding_positions['Mapping_position'] = 'Noncoding'
                    writefile(prot_id, out_dir, pident, isoform, consequence,
                              noncoding_positions, 'NoncodingPositions', csv, hdf)
                    left_positions = left_positions.loc[coding_positions_index]
                    del(noncoding_positions, coding_positions_index)
                
                left_positions['Protein_accession'] = prot_id

                # remove protein position since we know there are no matching positions
                # and will reduce the maximum number of combinations of rows
                #pdb.id ensembl.prot.id temp.chain int.chain 
                psdf = psdf.loc[:, ['PDB_code', 'Protein_accession',
                                    'PDB_chain', 'Protein_alignment_start',
                                     'Protein_alignment_end', 'PDB_alignment_start',
                                      'PDB_alignment_end', 'Pident']]
                #psdf = psdf.iloc[:, np.r_[0:3, 4:9]]
                psdf.drop_duplicates(inplace=True)
                # merge rest of positions with protein structures
                unmapped_positions = left_positions.join(
                    psdf.set_index('Protein_accession'), on='Protein_accession', how='inner')
                # convert col to numeric to make comparison
                #unmapped_positions['Protein_position'] = pd.to_numeric(
                #    unmapped_positions['Protein_position'], errors='coerce')
                unmapped_positions['Protein_position'] = unmapped_positions.Protein_position.values.astype(
                    int)
                unmapped_positions[(unmapped_positions.Protein_position.values >= unmapped_positions.Protein_alignment_start.values.astype(int)) & (
                    (unmapped_positions.Protein_position <= unmapped_positions.Protein_alignment_end.astype(int)))]
                
                if unmapped_positions.empty is False:
                    cs = annovars.columns.values.tolist()
                    cs.append("Protein_accession")  
                    unmapped_positions = unmapped_positions[[c for c in unmapped_positions.columns if c in cs]]
                    unmapped_positions.drop_duplicates(inplace=True)
                    unmapped_positions['Mapping_position'] = 'Unmapped'
                    writefile(prot_id, out_dir, pident, isoform, consequence,
                            unmapped_positions, 'UnmappedPositions', csv, hdf)
                
                    del(left_positions, unmapped_positions)  
                # mapped position is on the rest of the structure
            structure_positions = mapped_positions[mapped_positions['Interaction_type'].isna()]
            structure_positions = structure_positions.drop(['Chimera_interacting_position', 'Chimera_3D_position',
                                                              'PDB_interacting_3D_position','PDB_interacting_aa',
                                                              'Interface_min_distance', 'PDB_interacting_B_factor',
                                                              'PDB_interacting_chain', 'Interaction_type'], axis=1, errors = 'ignore')
            #do proper arragenments if no resulst are retrieved
            if structure_positions.empty is False:
                structure_positions.drop_duplicates(inplace=True)
                #structure_positions['PDB_seq_position'] = structure_positions['Protein_position'] - \
                # structure_positions['Protein_alignment_start'] + structure_positions['PDB_alignment_start']
                structure_positions['Mapping_position'] = 'Structure'
                writefile(prot_id, out_dir, pident, isoform, consequence,
                        structure_positions, 'StructurePositions', csv, hdf)
                    #unmapped_positions = left_positions#[~left_positions.index.isin(structure_positions.index)]

                # if unmapped_positions.empty is False:
                #     cs = annovars.columns.values.tolist()
                #     cs.append("Protein_accession")  
                #     unmapped_positions = unmapped_positions[[c for c in unmapped_positions.columns if c in cs]]
                #     unmapped_positions.drop_duplicates(inplace=True)
                #     unmapped_positions['Mapping_position'] = 'Unmapped'
                #     writefile(prot_id, out_dir, pident, isoform, consequence,
                #             unmapped_positions, 'UnmappedPositions', csv, hdf)
                # else:
                #     try:
                #         unmapped_positions = left_positions
                #         if unmapped_positions.empty is False:
                #             unmapped_positions.drop_duplicates(inplace=True)
                #             unmapped_positions['Mapping_position'] = 'Unmapped'
                #             writefile(prot_id, out_dir, pident, isoform, consequence,
                #                       unmapped_positions, 'Unmappedpositions', csv, hdf)
                #     except:
                #         pass
                del (structure_positions)
        ###########################################################################
        # stop if there are no results
        if mapped_positions.empty:
            # report results
            logger.warning('Warning: ' + prot_id +
                           ' does not map with any annotated position.\n')
        # if merging was successful, create setID file and
        # save the merged dataframe as well
        else:
            mapped_positions['Mapping_position'] = 'Interface'
            mapped_positions = mapped_positions[mapped_positions['Interaction_type'].notna()]
            setID_file = mapped_positions[['Structure_feature_id',
                                          'Uploaded_variation']]
            setID_file.drop_duplicates(inplace=True)
            mapped_positions.drop_duplicates(inplace=True)
            # Save the merged dataframe, appending results and not
            #  reapeting headers
            with open(os.path.join(out_dir, ('setID_pident' + str(pident) + '_isoform_' +
                                             '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '.txt')), 'a') as f:
                setID_file.to_csv(f, sep=',', index=False,
                                  header=f.tell() == 0)
            writefile(prot_id, out_dir, pident, isoform, consequence,
                      mapped_positions, 'InterfacePositions', csv, hdf)
            del(mapped_positions)
        del(psdf, annovars)
