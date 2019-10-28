#!/usr/bin/env python3
import sys
import os
import re
import pandas as pd
import numpy as np

def translate_ensembl(biomartdb, ensid):
    
    with open(biomartdb) as f:

        cols = f.readline()
      
        # find all the lines in the VEP file that contain information for a
        #  certain geneID
        for line in f:

            if ensid in line:  # it is faster than regex
                
                gene_col = cols.index("Gene stable ID")
                prot_col = cols.index("Protein stable ID\n")
                if "ENSP" in ensid:
                    protID = ensid
                    # remove \n at the end of the word
                    geneID = line.split(",")[gene_col].strip()
                    return {'protID': protID, 'geneID': geneID}
                elif "ENSG" in ensid:
                    protID = line.split(",")[prot_col].strip()
                    geneID = ensid
                    return {'protID': protID, 'geneID': geneID}
                else:
                    return IOError('Input not recognized')
                    #return
       
        return IOError('Input not recognized')


        