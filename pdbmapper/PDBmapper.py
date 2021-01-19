# coding: utf-8

# Import necesary modules
import sys
import os
import re
import glob
import pandas as pd
import numpy as np
from .db_parser import parser
# .interface_parser import reshape
from .decorator import tags
from .explode import explode
from .logger import get_logger


def PDBmapper(protid,  geneid, transcritpID, psdb, vardb, out_dir, pident, isoform, APPRIS, UniprotID, consequence, loc, varid=None, hdf = False):
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

    # parse interfaces corresponding to the selected protein ID
    psdf = parser(protid, psdb)

    if psdf.empty:
        logger.error('Interfaces of protein ' +
                     protid + 'could not be parsed.')
        raise IOError
    else:
        logger.info('Interfaces of protein ' + protid + ' parsed.')

    # if database is compacted and explode is needed

    columns_type = psdf.applymap(lambda x: isinstance(x, list)).all()
    columns_list = columns_type.index[columns_type].tolist()

    # cols_stack
    cols_stack = psdf.apply(lambda x: x.astype(
        str).str.match(r'[a-zA-Z0.-9]+-[a-zA-Z0.-9]+'))
    colsnames_stack = psdf.columns[cols_stack.any()].tolist()

    if any(cols_stack) is True:
        for col in colsnames_stack:

            psdf.loc[:, col] = psdf[col].str.replace(
                '-', ',')
            psdf.loc[:, col] = psdf[col].str.split(',')

        psdf = explode(
            psdf, colsnames_stack)

    elif columns_list is not None:
        psdf = explode(
            psdf, columns_list)
    else:
        psdf[colsnames_stack] = \
            psdf[colsnames_stack].astype(str)

    # if default database is used minor modifications are needed
    if pident is not None:
        logger.info('Filtering interfaces by pident = ' + str(pident) + '%.')
        # filter by pident
        pident = int(pident)  # from str to int
        psdf_pident = psdf.loc[psdf.Pident >= pident]
        # if pident threshold is to high, the next maximum value of pident is
        # notified in log file
        if psdf_pident.empty:
            alt_pident = psdf.loc[:, "Pident"].max()
            logger.error('Warning: for protid ' + str(pident) +
                         ', the variable "Pident" equal to ' +
                         str(pident) + ' is too high.\n A threshold lower than or equal to ' +
                         str(alt_pident) + ' would retrieve results.')

            raise IOError()

    # parse variants corresponding to the selected protein ID
    annovars = parser(geneid, vardb)
    logger.info('Variants file from gene id ' + geneid + ' parsed.')

    # filter by transcript ID
    annovars = annovars[annovars['Feature'] == transcritpID]
    if annovars.empty:
        raise IOError()

    # filter by variant type if one or more selected
    if consequence is not None:
        annovars = annovars[annovars['Consequence'].astype(
            str).str.contains('|'.join(consequence))]
        logger.info('Filter of features = ' + str(consequence))

        # if filter returns an empty df, raise error
        if annovars.empty is True:
            logger.error(
                'Variants could not be filtered by feature type = ' + consequence)
            raise IOError()

    # filter by variant type if one or more selected
    if varid is not None:
        if 'Existing_variation' in annovars.columns:
            annovars = annovars[
                annovars['Uploaded_variation'].astype(
                    str).str.contains(varid) |
                annovars['Existing_variation'].astype(
                    str).str.contains(varid)]
        else:
            annovars = annovars[
                annovars['Uploaded_variation'].astype(
                    str).str.contains(varid)]
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
        sub_df = explode(sub_df, ['Protein_position'])
        # drop unnecesary columns
        sub_df.drop(['start', 'end'], inplace=True, axis=1)
        # concatenate final result
        annovars = pd.concat([remaining_df, sub_df], sort=False)

    # for sucessful merge, Protein_position column must be str type
    psdf['Protein_position'] = psdf['Protein_position'].astype(str)
    annovars['Protein_position'] = annovars['Protein_position'].astype(str)
    annovars['APPRIS_isoform'] = APPRIS
    annovars['Uniprot_accession'] = UniprotID
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
            left_variants = annovars.drop(mapped_variants.index)
        except:
            left_variants = annovars
        # remove non protein coding variants
        left_variants.drop_duplicates(inplace=True)
        noncoding_variants_index = left_variants.Protein_position.str.contains(
            '\.|\-', regex=True)

        noncoding_variants = left_variants.loc[noncoding_variants_index]

        # non-protein coding mutations
        if noncoding_variants.empty is False:
            noncoding_variants = noncoding_variants.drop(
                columns=['Uniprot_accession'])  # not needed
            noncoding_variants['Mapping_position'] = 'Noncoding'
            with open(os.path.join(out_dir, ('NoncodingVariants_pident' + str(pident) + '_isoform_' +
                                             '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '.File')), 'a') as f:
                noncoding_variants.to_csv(f, sep=',', index=False,
                                          header=f.tell() == 0)

        left_variants = left_variants.loc[~noncoding_variants_index]
        left_variants['Protein_accession'] = protid

        # remove protein position since we know there are no matching positions
        # and will reduce the maximum number of combinations of rows
        psdf = psdf.iloc[:, np.r_[0:3, 4:9]]
        psdf.drop_duplicates(inplace=True)
        # merge rest of variants with protein structures
        left_variants = left_variants.merge(
            psdf.set_index('Protein_accession'), on='Protein_accession', how='inner')
        # convert col to numeric to make comparison
        left_variants['Protein_position'] = pd.to_numeric(
            left_variants['Protein_position'], errors='coerce')
        left_variants['Protein_position'] = left_variants.Protein_position.values.astype(
            int)
        # mapped variant is on the rest of the structure
        structure_variants = left_variants[(left_variants.Protein_position.values >= left_variants.Protein_start_position.values) & (
            (left_variants.Protein_position <= left_variants.Protein_end_position))]

        # do proper arragenments if no resulst are retrieved
        if structure_variants.empty is False:

            if mapped_variants.empty is False:
                mapped_variants['Mapping_position'] = 'Interface'

            structure_variants.drop_duplicates(inplace=True)
            structure_variants['Mapping_position'] = 'Structure'
            with open(os.path.join(out_dir, ('StructureVariants_pident' + str(pident) + '_isoform_' +
                                             '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '.File')), 'a') as f:
                structure_variants.to_csv(f, sep=',', index=False,
                                          header=f.tell() == 0)
            # unmapped variants
            unmapped_variants = left_variants.drop(structure_variants.index)
            if unmapped_variants.empty is False:
                unmapped_variants = unmapped_variants.iloc[:, 0:16]
                unmapped_variants.drop_duplicates(inplace=True)
                unmapped_variants['Mapping_position'] = 'Unmapped'
                with open(os.path.join(out_dir, ('UnmappedVariants_pident' + str(pident) + '_isoform_' +
                                                 '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '.File')), 'a') as f:
                    unmapped_variants.to_csv(f, sep=',', index=False,
                                             header=f.tell() == 0)
        else:
            if mapped_variants.empty is False:
                mapped_variants['Mapping_position'] = 'Interface'
                unmapped_variants = left_variants.drop(mapped_variants.index)

                if unmapped_variants.empty is False:
                    unmapped_variants = unmapped_variants.iloc[:, 0:16]
                    unmapped_variants.drop_duplicates(inplace=True)
                    unmapped_variants['Mapping_position'] = 'Unmapped'
                    with open(os.path.join(out_dir, ('UnmappedVariants_pident' + str(pident) + '_isoform_' +
                                                     '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '.File')), 'a') as f:
                        unmapped_variants.to_csv(f, sep=',', index=False,
                                                 header=f.tell() == 0)
            else:
                # report results
                logger.warning('Warning: ' + protid +
                               ' does not map with any annotated variant.\n')
                # unmapped variants
                unmapped_variants = left_variants
                if unmapped_variants.empty is False:
                    unmapped_variants = unmapped_variants.iloc[:, 0:16]
                    unmapped_variants.drop_duplicates(inplace=True)
                    unmapped_variants['Mapping_position'] = 'Unmapped'
                    with open(os.path.join(out_dir, ('UnmappedVariants_pident' + str(pident) + '_isoform_' +
                                                     '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '.File')), 'a') as f:
                        unmapped_variants.to_csv(f, sep=',', index=False,
                                                 header=f.tell() == 0)
                raise IOError()
        del (structure_variants,
             unmapped_variants, left_variants, noncoding_variants, noncoding_variants_index)
    ###########################################################################

    # stop if there are no results
    if mapped_variants.empty:
        # report results
        logger.warning('Warning: ' + protid +
                       ' does not map with any annotated variant.\n')
        raise IOError()

    # if merging was successful, create setID file and
    # save the merged dataframe as well
    else:
        setID_file = mapped_variants[['Structure_feature_id',
                                      'Uploaded_variation']]
        setID_file.drop_duplicates(inplace=True)
        mapped_variants.drop_duplicates(inplace=True)

        # Save the merged dataframe, appending results and not
        #  reapeting headers
        with open(os.path.join(out_dir, ('setID_pident' + str(pident) + '_isoform_' +
                                         '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '.File')), 'a') as f:
            setID_file.to_csv(f, sep=',', index=False,
                              header=f.tell() == 0)

        with open(os.path.join(out_dir, ('MappedVariants_pident' + str(pident) + '_isoform_' +
                                         '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '.File')), 'a') as f:

            mapped_variants.to_csv(f, sep=',', index=False,
                                   header=f.tell() == 0)
        if hdf is True: 
            fn = os.path.join(out_dir, ('MappedVariants_pident' + str(pident) + '_isoform_' +
                                         '_'.join(isoform) + '_consequence_' + '_'.join(consequence) + '.File.hdf5'))
            store = pd.HDFStore(fn)

            store.append('mapped_variants', mapped_variants, format='t',  data_columns=True)
            

            #mapped_variants.to_hdf(f, key = 'mapped_variants', mode = 'r+', append = True, format = 't')

    del(psdf, psdf_pident, annovars)
