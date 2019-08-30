#!/usr/bin/env python3

#import sys
#import os
#import re
#import pandas as pd
#from timeit import default_timer as timer

# from vcf_parser import (Genotype, HeaderParser)
# from vcf_parser.utils import (format_variant, split_variants)


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
    print(vepf)


    
def VEP_parser( VEPfile, geneID):
    # create empty list to store the rows 
    matches = []
    for line in VEPfile:
        df = re.findall(r''+geneID, line) # similar to grep. Faster than reading the 
        if df: 
            matches.append(line.split("\t"))
    df0 = pd.read_csv(VEPfile, nrows = 0, skiprows = 42, sep = "\t")
    df = pd.DataFrame(matches)
    df.columns = df0.columns

    #index = df_vep.columns.get_loc("Amino_acids")
    return df


def ProtIntDB_parser(InterfacesDB, protID):
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
    colnames = InterfacesDB.readline().split(" ")
    df.columns = colnames
    
    return df

