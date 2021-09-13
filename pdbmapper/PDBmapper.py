# coding: utf-8

# Import necesary modules
import sys
import os
import re
import glob
import pandas as pd
import numpy as np
import dask.dataframe as dd

from .db_parser import parser
# .interface_parser import reshape
from .decorator import tags
from .explode import explode
from .explode2 import explode2
from .logger import get_logger
from .writefile import writefile

#import vaex


def PDBmapper(protid,  geneid, transcriptid, psdb, vardb, out_dir, pident, evalue, isoform, APPRIS, consequence, loc, varid=None, csv=False, hdf=False):
    '''
    Map interfaces and genomic anntoated variants and returns a
    setID.File, necessary input for SKAT. Additionaly, it creates
    another file with detailed information regarding the maped areas.

    Parameters
    ----------
    protid : str
        Ensembl protein ID
    geneid : str
        Translated Ensembl protein ID
    psdb : str
        Directory where to find interface database
    vardb : str
        Directory where to find variants database
    out_dir : str
        Output directory
    pident : int
        Thershold of sequence identity (percertage).

    Returns
    -------
    setID.File
        txt file containing a data frame two columns corresponding to the
        analyzed interface id and the corresponding annotated genomic variants.
    MappedVariants.File
        Same as setID.File but with additional information describing the
        interfaces and the variants.
    '''
    # log file
    logger = get_logger(' PDBmapper', out_dir)
    # parse variants corresponding to the selected protein ID
    try:
        annovars = parser(transcriptid, vardb)     
        if consequence is not None:
            annovars = annovars[annovars['Consequence'].astype(
                str).str.contains('|'.join(consequence))]
            logger.info('Filter of features = ' + str(consequence))

            # if filter returns an empty df, raise error
            if annovars.empty is True:
                logger.error(
                    'Variants could not be filtered by feature type = ' + consequence)
                raise IOError()
        else:
            consequence = ['all']

        if isoform is None:
            isoform = ['all']

        # filter by variant type if one or more selected
        if varid is not None:
            if 'Existing_variation' in annovars.columns:
                annovars = annovars[
                    (annovars.Uploaded_variation.astype(str) == str(varid)) |
                    (annovars.Existing_variation.astype(str) == str(varid))]
            else:
                annovars = annovars[
                    (annovars.Uploaded_variation.astype(str) == str(varid))]
            logger.info('Variant \'' + str(varid) + '\' has been selected.')
            # if filter returns an empty df, raise error
            if annovars.empty:
                logger.error(
                    'Variants could not be filtered by variant id \'' + str(varid) + '\'')
                raise IOError()
        # for variants with high impact affecting several aminoacidic positions,
        # the protein position is a range. split the range to have each position
        # individually
        if any(annovars['Protein_position'].astype(str).str.contains(r'[0-9]-[0-9]')):
            # subset hight impact variants
            sub_df = annovars[annovars['Protein_position'].astype(str).str.contains(
                r'[0-9]-[0-9]')]
            # subset the remaining variants to concatenate afterwards
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
        psdf = parser(protid, psdb)
        logger.info('Protein features file of ' + protid + ' parsed.')
        # Get col PDB_position
        # if database is compacted and explode is needed
        columns_type = psdf.applymap(lambda x: isinstance(x, list)).all()
        columns_list = columns_type.index[columns_type].tolist()

        # cols_stack
        cols_stack = psdf.apply(lambda x: x.astype(
            str).str.match(r'[a-zA-Z0.-9]+-[a-zA-Z0.-9]+'))
        colsnames_stack = psdf.columns[cols_stack.any()].tolist()
        # add column for chimera script
        
        if 'PDB_interacting_3D_position' in colsnames_stack:
            #psdf['Interface_interacting_positions'] = psdf['PDB_interacting_position']  
            psdf['Chimera_interacting_position'] = psdf['PDB_interacting_3D_position']
        if 'Evalue' in colsnames_stack:
            colsnames_stack.remove('Evalue')
        if 'PDB_code' in colsnames_stack:
            colsnames_stack.remove('PDB_code')
        if 'PDB_3D_position' in colsnames_stack:
            #colsnames_stack.remove('PDB_3D_position')
            psdf['Chimera_3D_position'] = psdf['PDB_3D_position']
        if 'Structure_feature_id' in colsnames_stack:
            colsnames_stack.remove('Structure_feature_id')

        if any(colsnames_stack):
            psdf = explode(psdf , colsnames_stack, '-')
        elif any(columns_list):
            psdf = explode(
                psdf, columns_list, '-')
        else:
            psdf[colsnames_stack] = \
                psdf[colsnames_stack].astype(str)
 
        # if default database is used minor modifications are needed
        if pident is not None:
            logger.info('Filtering interfaces by pident = ' +
                        str(pident) + '%.')
            # filter by pident
            pident = int(pident)  # from str to int
            psdf = psdf.loc[psdf.Pident >= pident]
            # if pident threshold is to high, the next maximum value of pident is
            # notified in log file
            if psdf.empty:
                alt_pident = psdf.loc[:, "Pident"].max()
                logger.error('Warning: for protid ' + str(pident) +
                             ', the variable "Pident" equal to ' +
                             str(pident) + ' is too high.\n A threshold lower than or equal to ' +
                             str(alt_pident) + ' would retrieve results.')

                raise IOError()
        if evalue is not None:
            logger.info('Filtering interfaces by evalue = ' +
                        str(evalue) + '%.')
            # filter by pident
            evalue = float(evalue)  # from str to int
            psdf = psdf.loc[psdf.Evalue >= evalue]
            # if pident threshold is to high, the next maximum value of pident is
            # notified in log file
            if psdf.empty:
                alt_evalue = psdf.loc[:, "Evalue"].min()
                logger.error('Warning: for protid ' + str(evalue) +
                             ', the variable "Evalue" equal to ' +
                             str(evalue) + ' is too low.\n A threshold higher than or equal to ' +
                             str(alt_evalue) + ' would retrieve results.')

                raise IOError()       
    except IOError:
        psdf = False
        
    if psdf is False and annovars is False:
        logger.error('Protein ' +
                     protid + 'could not be parsed.')
        raise IOError

    elif psdf is not False and annovars is False:
        logger.error('Protein ' +
                     protid + 'has no mapping variants.')
        raise IOError

    elif psdf is False and annovars is not False:

        if loc:
            # remove non protein coding variants
            annovars.drop_duplicates(inplace=True)
            try:
                noncoding_variants_index = annovars.Amino_acids.str.contains(
                    '\.|\-', regex=True, na=False)
                noncoding_variants = annovars.loc[noncoding_variants_index]
            except:
                noncoding_variants = False
            # non-protein coding mutations
            if noncoding_variants is not False:
                noncoding_variants['APPRIS_isoform'] = ''
                noncoding_variants['Mapping_position'] = 'Noncoding'
                writefile(protid, out_dir, pident, isoform, consequence,
                          noncoding_variants, 'NoncodingVariants', csv, hdf)
                unmapped_variants = annovars.loc[~noncoding_variants_index]
            else:
                unmapped_variants = annovars

            if unmapped_variants.empty is False:
                unmapped_variants.drop_duplicates(inplace=True)
                unmapped_variants['APPRIS_isoform'] = ''
                unmapped_variants['Mapping_position'] = 'Unmapped'
                writefile(protid, out_dir, pident, isoform, consequence,
                          unmapped_variants, 'UnmappedVariants', csv, hdf)
        
    elif psdf is not False and annovars is not False:
        # for sucessful merge, Protein_position column must be str type
        psdf['Protein_position'] = psdf['Protein_position'].astype(str)
        annovars['Protein_position'] = annovars['Protein_position'].astype(str)
        annovars['APPRIS_isoform'] = APPRIS
       # annovars['Uniprot_accession'] = UniprotID
        # Merge them both files
        mapped_variants = annovars.join(
            psdf.set_index('Protein_position'), on='Protein_position', how='inner')

        if isoform is None:
            isoform = ['all']
        if consequence is None:
            consequence = ['all']
        ###########################################################################
        # Locate rest of variants (mapping to a structure or not)
        ###########################################################################
        if loc:
            # remove already mapped variants
            try:
                left_variants = annovars.drop(list(set(mapped_variants.index)))
            except:
                left_variants = annovars
            # remove non protein coding variants
            if left_variants.empty is False:
                left_variants.drop_duplicates(inplace=True)
                noncoding_variants_index = left_variants.Amino_acids.str.contains(
                    '\.|\-', regex=True, na=False)

                noncoding_variants = left_variants.loc[noncoding_variants_index]

                # non-protein coding mutations
                if noncoding_variants.empty is False:
                   # noncoding_variants = noncoding_variants.drop(
                   #     columns=['Uniprot_accession'])  # not needed
                    noncoding_variants['Mapping_position'] = 'Noncoding'
                    writefile(protid, out_dir, pident, isoform, consequence,
                              noncoding_variants, 'NoncodingVariants', csv, hdf)
                    left_variants = left_variants.loc[~noncoding_variants_index]

                left_variants['Protein_accession'] = protid
                # remove protein position since we know there are no matching positions
                # and will reduce the maximum number of combinations of rows
                #pdb.id ensembl.prot.id temp.chain int.chain 
                psdf = psdf.loc[:, ['PDB_code', 'Protein_accession',
                                    'PDB_chain', 'Protein_alignment_start',
                                     'Protein_alignment_end', 'PDB_alignment_start',
                                      'PDB_alignment_end', 'Pident']]
                #psdf = psdf.iloc[:, np.r_[0:3, 4:9]]
                psdf.drop_duplicates(inplace=True)
                # merge rest of variants with protein structures
                left_variants = left_variants.join(
                    psdf.set_index('Protein_accession'), on='Protein_accession', how='inner')
                # convert col to numeric to make comparison
                left_variants['Protein_position'] = pd.to_numeric(
                    left_variants['Protein_position'], errors='coerce')
                left_variants['Protein_position'] = left_variants.Protein_position.values.astype(
                    int)
                # mapped variant is on the rest of the structure
                structure_variants = left_variants[(left_variants.Protein_position.values >= left_variants.Protein_alignment_start.values) & (
                    (left_variants.Protein_position <= left_variants.Protein_alignment_end))]

                # do proper arragenments if no resulst are retrieved
                if structure_variants.empty is False:
                    structure_variants.drop_duplicates(inplace=True)
                    
                    # structure_variants['PDB_position'] = np.where(structure_variants['Protein_alignment_start'] > structure_variants['PDB_alignment_start'],
                    #  structure_variants['Protein_position'] - structure_variants['Protein_alignment_start'] + structure_variants['PDB_alignment_start'],
                    #   np.where(structure_variants['Protein_alignment_start'] < structure_variants['PDB_alignment_start'],
                    #   structure_variants['Protein_position'] + structure_variants['PDB_alignment_start'] -1,
                    #   structure_variants['Protein_position']  ))

                    structure_variants['PDB_seq_position'] = structure_variants['Protein_position'] - structure_variants['Protein_alignment_start'] + structure_variants['PDB_alignment_start']

                    structure_variants['Mapping_position'] = 'Structure'
                    writefile(protid, out_dir, pident, isoform, consequence,
                              structure_variants, 'StructureVariants', csv, hdf)
                    # unmapped variants
                   # unmapped_variants = left_variants.drop(
                   #     structure_variants.index)
                    unmapped_variants = left_variants[~left_variants.index.isin(structure_variants.index)]
                   # unmapped_variants = left_variants.loc[left_variants.index.drop(structure_variants.index)]

                    if unmapped_variants.empty is False:
                        cs = annovars.columns.values.tolist()
                        cs.append("Protein_accession")  
                        unmapped_variants = unmapped_variants[[c for c in unmapped_variants.columns if c in cs]]
                        unmapped_variants.drop_duplicates(inplace=True)
                        unmapped_variants['Mapping_position'] = 'Unmapped'
                        writefile(protid, out_dir, pident, isoform, consequence,
                                  unmapped_variants, 'UnmappedVariants', csv, hdf)

                else:
                    try:
                        unmapped_variants = left_variants
                        if unmapped_variants.empty is False:
                            unmapped_variants.drop_duplicates(inplace=True)
                            unmapped_variants['Mapping_position'] = 'Unmapped'
                            writefile(protid, out_dir, pident, isoform, consequence,
                                      unmapped_variants, 'UnmappedVariants', csv, hdf)
                    except:
                        pass

                del (structure_variants, left_variants,
                     noncoding_variants, noncoding_variants_index)  
        ###########################################################################
        # stop if there are no results

        if mapped_variants.empty:
            # report results
            logger.warning('Warning: ' + protid +
                           ' does not map with any annotated variant.\n')
            #raise IOError()

        # if merging was successful, create setID file and
        # save the merged dataframe as well
        else:
            
            mapped_variants['Mapping_position'] = 'Mapped'
            setID_file = mapped_variants[['Structure_feature_id',
                                          'Uploaded_variation']]
            setID_file.drop_duplicates(inplace=True)
            mapped_variants.drop_duplicates(inplace=True)

            # Save the merged dataframe, appending results and not
            #  reapeting headers
            with open(os.path.join(out_dir, ('setID_pident' + str(pident) + '_isoform_' +
                                             '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '.txt')), 'a') as f:
                setID_file.to_csv(f, sep=',', index=False,
                                  header=f.tell() == 0)
            writefile(protid, out_dir, pident, isoform, consequence,
                      mapped_variants, 'MappedVariants', csv, hdf)
            # or 'batch_1.arrow'

            #store = pd.HDFStore(fn)
            #store.append('mapped_variants', mapped_variants, format='t',  data_columns=True)
            #mapped_variants.to_hdf(f, key = 'mapped_variants', mode = 'r+', append = True, format = 't')
            del(mapped_variants)
        del(psdf, annovars)
