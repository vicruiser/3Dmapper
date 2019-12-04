# -*- coding: utf-8 -*-

# import necessary modules
import sys
import os
import re
import pandas as pd


def translate_ensembl(ensid, isoform_filter):
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
    dirname = os.path.dirname(__file__)
    # rel_path =
    biomartdb = os.path.join(dirname, "data/biomart_GRCh38p13_nov2019.dat")

    with open(biomartdb) as f:
        # get col names
        cols = f.readline().split(',')
        # parse file
        for line in f:
            # stop when find input ensembl id
            if ensid in line:  # it is faster than regex
                # detect wether the input is protein or gene id.
                gene_col = cols.index("Gene stable ID")
                prot_col = cols.index("Protein stable ID\n")
                transcript_col = cols.index("Transcript stable ID")
                # check type of isoform
                isoform_col = cols.index("APPRIS annotation")
                isoform_class = line.split(",")[isoform_col].strip()
                # print(isoform_col)
                if isoform_class in isoform_filter:
                    if "ENSP" in ensid:
                        protID = ensid
                        geneID = line.split(",")[gene_col].strip()
                        transcriptID = line.split(",")[transcript_col].strip()
                        return {'protID': protID, 'geneID': geneID, 'transcriptID': transcriptID}
                    elif "ENSG" in ensid:
                        protID = line.split(",")[prot_col].strip()
                        geneID = ensid
                        transcriptID = line.split(",")[transcript_col].strip()
                        return {'protID': protID, 'geneID': geneID, 'transcriptID': transcriptID}
                    # error handling
                    else:
                        raise IOError()
                else:
                    raise IOError()
        # error handling
        raise IOError()
