#!/usr/bin/env python3
#from __future__ import print_function

import sys
import os
import gzip
import re
import pandas as pd
#import dask.dataframe as dd
#import re
#import pkg_resources
#import click
#import locale

# from codecs import open, getreader


# from vcf_parser import (Genotype, HeaderParser)
# from vcf_parser.utils import (format_variant, split_variants)

####            Parser:         ####

# class VEPparser(object):

# def joinDataFrames(vcf_file, interfaces_db):
#     # Read the CSVs
#     df1 = dd.read_csv(vcf_file)
#     df2 = dd.read_csv(interfaces_db)

#     # Merge them
#     df = dd.merge(df1, df2, on='Bin_ID').compute()

#     # Save the merged dataframe
#     df.to_csv('merged.csv', index=False)

VEPfile = open("/home/vruizser/PhD/2018-2019/PanCancer_data/vep_output/vep_PanCan_chr_1_1-100000", "r")
geneID = 'ENSG00000227232'

def VEPparser(VEPfile, geneID):
    # create empty list to store 
    matches = []
    for line in VEPfile:
        df = re.findall(r''+geneID, line)
        if df: 
            matches.append(line.split("\t"))

    df = pd.DataFrame(matches)
    
    print(df[10]) 

def InterfacesDBparser(InterfacesDB, geneID):
    # create empty list to store 
    matches = []
    # find all the lines in the VEP file that contain information for a certain geneID
    for line in VEPfile:
        df = re.findall(r''+geneID, line)
        #if we find a match, store it
        if df:  
            matches.append(line.split("\t"))

    df = pd.DataFrame(matches)
    
    print(df[10]) 
    return

# define main function to execute the previous defined functions together
def main():

    VEPparser(VEPfile, geneID)
    return


##########################
# execute main function  #
##########################
if __name__ == '__main__':
    main()