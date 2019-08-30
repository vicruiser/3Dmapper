#!/usr/bin/env python3
#from __future__ import print_function

#import sys
#import os
#import gzip
#import re
#import pandas as pd
#from parse import Lexer, Parser, Token, State, NFA, Handler
#import dask.dataframe as dd
#import re
#import pkg_resources
#import click
#import locale
import VEPparser as vp
import sys
import os
import gzip
import re
import pandas as pd
from timeit import default_timer as timer
# from codecs import open, getreader


# from vcf_parser import (Genotype, HeaderParser)
# from vcf_parser.utils import (format_variant, split_variants)

####            Parser:         ####


def ensemblIDtranslator(biomartdb, ensid): 
    cols = biomartdb.readline().split(",")
    # find all the lines in the VEP file that contain information for a certain geneID
    for line in biomartdb:

        if ensid in line:  # it is faster than regex
            
            gene_col = cols.index("Gene stable ID")
            prot_col = cols.index("Protein stable ID\n")
            
            if "ENSP" in ensid:
                protID = ensid
                geneID = line.split(",")[gene_col].strip() # remove \n at the end of the word
            elif "ENSG" in ensid:
                protID = line.split(",")[prot_col].strip()
                geneID = ensid
            else :
                print("wrong id format")

    #df = pd.DataFrame(matches)
    return {'protID': protID, 'geneID': geneID}

def getVEPfile (crossref_file, geneID): 
     
    matches = []
    for line in crossref_file:
        if geneID in line: 
            vepf = line.split("\t")[1].strip()
    print(vepf)


def VEPparser( VEPfile, vepfile, geneID):
    # create empty list to store the rows 
    matches = []
    for line in VEPfile:
        df = re.findall(r''+geneID, line) # similar to grep. Faster than reading the 
        if df: 
            matches.append(line.split("\t"))
    df0 = pd.read_csv(vepfile, nrows = 0, skiprows = 42, sep = "\t")
    df = pd.DataFrame(matches)
    df.columns = df0.columns

    #index = df_vep.columns.get_loc("Amino_acids")
    return df


def InterfacesDBparser(InterfacesDB, protID):
    # create empty list to store 
    matches = []
    # find all the lines in the VEP file that contain information for a certain geneID
    for line in InterfacesDB:
        #df = re.findall(r'ENSP00000352835', line)
        if protID in line:  # it is faster than regex
            df = line
        #if we find a match, store it
            if df:  
                matches.append(line.split(" "))
    df = pd.DataFrame(matches)
    colnames = InterfacesDB.readline()
    df.columns = colnames
    
    return df


# def joinDataFrames(vcf_file, interfaces_db):
#     # Read the CSVs
#     df1 = dd.read_csv(vcf_file)
#     df2 = dd.read_csv(interfaces_db)

#     # Merge them
#     df = dd.merge(df1, df2, on='Bin_ID').compute()

#     # Save the merged dataframe
#     df.to_csv('merged.csv', index=False)
