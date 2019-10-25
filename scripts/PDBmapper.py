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

@tags(text_start = "Runing PDBmapper...",
      text_succeed = "Runing PDBmapper...done.",
      text_fail = "Runing PDBmapper...failed!",
      emoji = "🦸")
def PDBmapper(pid, geneid, int_db_dir, vcf_db_dir, out_dir):
    '''Generate setID file.
  
    Parameters
    ----------
    protID : str
        input path of the VEP file
    interfacesDB : str
        chosen name of the output file
    VCF_subset : str
        chosen name of the output file

    Returns
    -------
    setID.File
        write data frame to a txt file with two columns. One is the gene ids \
             and the other one the VEP file
    MappedVariants.File
        more into
    '''
    
    annovars = parser(geneid, vcf_db_dir, "\t")
    annoint  = parser(pid, int_db_dir, " ")
    annoint_reshape = reshape(annoint)

    # Merge them both files
    mapped_variants = pd.concat([annovars, annoint_reshape],
                                axis=1, join='inner')
    # stop if there are no results
    if mapped_variants.empty:
        print('Warning:', pid, 'does not map with any annotated variant.')
        exit(-1)
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


