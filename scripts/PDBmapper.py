#!/usr/bin/env python3
# coding: utf-8

# Import necesary modules
#import mapping_tools as mt
import numpy as np
import sys
import os
import gzip
import re
import pandas as pd
import glob
from scripts.db_parser import parser
from scripts.interface_parser import reshape
from scripts.decorator import tags


def PDBmapper(protID, geneID, int_db_dir, vcf_db_dir, out_dir, pident):
    '''
    Map interfaces and genomic anntoated variants and returns a
    setID.File, necessary input for SKAT. Additionaly, it creates
    another file with detailed information regarding the maped areas. 
  
    Parameters
    ----------
    protID : str
        Ensembl protein ID
    geneID : str
        Translated Ensembl protein ID
    int_db_dir : str
        Directory where to find interface database
    vcf_db_dir : str
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
    # parse interfaces corresponding to the selected protein ID
    annoint = parser(protID, int_db_dir, " ")
    # filter by pident
    pident = int(pident) # from str to int
    annoint_pident  = annoint.loc[annoint.pident >= pident]
    # if pident threshold is to high, the next maximum value of pident is set
    if annoint_pident.empty:
        # set alternative pident
        alt_pident = annoint.loc[: , "pident"].max()
        # register the change
        log = open(out_dir + '/log.File', 'a')
        log.write('Warning: for protID ' + protID +
                  ', the variable "pident" equal to ' +
                   pident + ' is too high.\n A threshold of ' + 
                   alt_pident + ' has been set to run PDBmapper, which is the maximum possible value.')
        # set the new subset
        annoint_pident  = annoint.loc[annoint.pident >= alt_pident]
    # spread the data frame to have one amino acid position per row instead of compacted.         
    annoint_reshape = reshape(annoint_pident)
    # parse variants corresponding to the selected protein ID
    annovars = parser(geneID, vcf_db_dir, " ")
    # Merge them both files
    mapped_variants = pd.concat([annovars, annoint_reshape],
                                axis=1, join='inner')
    # stop if there are no results
    if mapped_variants.empty:
        log = open(out_dir + '/log.File', 'a')
        log.write('Warning: ' + protID + ' does not map with any annotated variant.\n')
        #log.flush()
        raise IOError('Warning: ' + protID + ' does not map with any annotated variant.')
    
    # if merging was successful, create setID file and
    # save the merged dataframe as well
    else:
        setID_file = mapped_variants[['region_id',
                                       '#Uploaded_variation']]
        setID_file = setID_file.drop_duplicates()
      
        # Save the merged dataframe, appending results and not
        #  reapeting headers
        with open(out_dir + '/setID.File', 'a') as f:
            setID_file.to_csv(f, sep=' ', index=False,  header=f.tell() == 0)
        with open(out_dir + '/MappedVariants.File', 'a') as f:
            mapped_variants.to_csv(f, sep=' ', index=False,
                                   header=f.tell() == 0)


