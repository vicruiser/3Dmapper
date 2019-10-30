#!/usr/bin/env python3

#import necessary modules
import sys
import os
import re
import pandas as pd

def translate_ensembl(ensid):
    '''
    gene-protein Ensembl id translator. 
     
    Parameters
    ----------
    ensid : str
        Input Ensembl ID corresponding to a gene or protein. 
        
    Returns
    -------
    dict 
        Dictionary containing input ID and its corresponding tranlation.  
    '''
    # read reference file of ensembl ids
    biomartdb = '/home/vruizser/PhD/2018-2019/git/PDBmapper/default_input_data/gene_transcript_protein_ens_ids.txt'
    with open(biomartdb) as f:
        # get col names
        cols = f.readline()
        # parse file 
        for line in f:
            # stop when find input ensembl id
            if ensid in line:  # it is faster than regex
                # detect wether the input is protein or gene id. 
                gene_col = cols.index("Gene stable ID")
                prot_col = cols.index("Protein stable ID\n")
                if "ENSP" in ensid:
                    protID = ensid
                    geneID = line.split(",")[gene_col].strip()
                    return {'protID': protID, 'geneID': geneID}
                elif "ENSG" in ensid:
                    protID = line.split(",")[prot_col].strip()
                    geneID = ensid
                    return {'protID': protID, 'geneID': geneID}
                # error handling
                else:
                    raise IOError()
        # error handling
        raise IOError()


        