# -*- coding: utf-8 -*-
# import necessary modules
import os
import os.path
import time
import subprocess
import sys
import csv
import glob
import re
import pandas as pd
import numpy as np

from .decorator import tags
from .logger import get_logger
# define request function to avoid repeating code
#maf = "/home/vruizser/head_mc3.maf"
#maf = '/home/vruizser/PhD/2018-2019/Immunity_interfaces_analysis/raw_data/mc3.v0.2.8.PUBLIC.maf'

# add decorator to main function


@tags(text_start="Converting MAF to VEP...This might take up some time...\n",
      text_succeed="Converting MAF to VEP...done.\n",
      text_fail="Converting MAF to VEP...failed!\n",
      emoji="\U0001F504 ")
def maf2vep(input_file, out_dir, out_file, overwrite, log_dir):
    # log file
    logger = get_logger('maf2vep', log_dir)
    logger.info('Converting maf file to vep format.')
    # set the name of 'uploaded_variants'
    upvar_cols = ["Chromosome", "Start_Position", "Reference_Allele", "Allele"]
    # set the name of 'location'
    loc_cols = ["Chromosome", "Start_Position"]
    # set the name of  the rest of columns included in a VEP file
    rest_of_cols = ["Allele",
                    "Gene",
                    "Feature",
                    "Feature_type",
                    "Consequence",
                    "cDNA_position",
                    "CDS_position",
                    "Protein_position",
                    "Amino_acids",
                    "Codons",
                    "Existing_variation"]
    # set header
    vep_header = ['Uploaded_variation', 'Location'] + rest_of_cols
    # read the MAF file line by line to modify each of them in situ
    with open(input_file) as f:
        # get col names
        cols = f.readline().split('\t')
        # take the index position of upvar_cols in the file
        uploaded_variation_cols = np.where(np.isin(cols, upvar_cols))
        # take the index position of loc_cols in the file
        location_cols = np.where(np.isin(cols, loc_cols))
        # and the rest of columns
        rest_cols = np.where(np.isin(cols, rest_of_cols))
        # open new file to write results
        with open(os.path.join(out_file), 'w', newline='', encoding='utf8') as csvfile:
            # use csv module
            writer = csv.writer(csvfile, delimiter='\t')
            # first row is vep header.
            writer.writerow(vep_header)
            # parse file
            for line in f:
                # split row into list and select the elements of the corresponding
                # indexes, previously defined
                row_splitted = np.asarray(line.split('\t'))
                Uploaded_variation = row_splitted[uploaded_variation_cols]
                Location = row_splitted[location_cols]
                rest = row_splitted[rest_cols]
                # merge previous list objects conveniently according to VEP file standards
                to_save = ['_'.join(Uploaded_variation[:-1]) + '/' +
                           Uploaded_variation[-1]] + [':'.join(Location)] + rest.tolist()
                # save results
                writer.writerow(to_save)

    logger.info('Conversion done successfully.')
