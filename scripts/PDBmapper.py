#!/usr/bin/env python3
# coding: utf-8

# Import necesary modules
import mapping_tools as mt
import numpy as np
import sys
import os
import gzip
import re
import pandas as pd
from timeit import default_timer as timer
from VEPcrossref import VEPfileCrossrefGenerator as cr
import parse_argv
import glob


def PDBmapper(protID, interfacesDB, VCF_subset, output_dir):
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

    # Merge them both files
    mapped_variants = pd.concat([VCF_subset, interfacesDB],
                                axis=1, join='inner')
    # stop if there are no results
    if mapped_variants.empty:
        print('Warning:', protID, 'does not map with any annotated variant.')
        exit(-1)
    # if merging was successful, create setID file and
    # save the merged dataframe as well
    else:
        setID_file = mapped_variants[['region_id',
                                      '#Uploaded_variation']]
        setID_file = setID_file.drop_duplicates()
        # Save the merged dataframe, appending results and not
        #  reapeting headers
        with open(output_dir + 'setID.File', 'a') as f:
            setID_file.to_csv(f, sep=' ', index=False,  header=f.tell() == 0)
        with open(output_dir + 'MappedVariants.File', 'a') as f:
            mapped_variants.to_csv(f, sep=' ', index=False,
                                   header=f.tell() == 0)


