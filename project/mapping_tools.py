#!/usr/bin/env python3

import sys
import os
import re
import pandas as pd
import numpy as np
#from timeit import default_timer as timer

# from vcf_parser import (Genotype, HeaderParser)
# from vcf_parser.utils import (format_variant, split_variants)

######## credit to @piRSquared from stackoverflow ######################################
def explode(df, columns):
    idx = np.repeat(df.index, df[columns[0]].str.len())
    a = df.T.reindex_axis(columns).values
    concat = np.concatenate([np.concatenate(a[i]) for i in range(a.shape[0])])
    p = pd.DataFrame(concat.reshape(a.shape[0], -1).T, idx, columns)
    return pd.concat([df.drop(columns, axis=1), p], axis=1).reset_index(drop=True)   
#######################################################################################

def parser( input_file, ensemblID, colnames, sep):
    
    # create empty list to store the rows 
    matches = []
    for line in input_file:
        if ensemblID in line: # similar to grep. Faster than reading the 
            matches.append(line.split(sep))
    
    df = pd.DataFrame(matches)
    df.columns = colnames

    return df

def ensemblID_translator(biomartdb, ensid): 
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

    return {'protID': protID, 'geneID': geneID}



def VEP_getter (crossref_file, geneID): 
    matches = []
    for line in crossref_file:
        if geneID in line: 
            vepf = line.split("\t")[1].strip()
    
    return vepf
