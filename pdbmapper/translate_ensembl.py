# -*- coding: utf-8 -*-

# import necessary modules
import sys
import os
import re
import pandas as pd
from .run_subprocess import call_subprocess
from .logger import get_logger


def translate_ensembl(ensid, log_dir, isoform_filter=None):
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
    # set logger
    logger = get_logger('translate_ensembl', log_dir)
    # read reference file of ensembl ids
    dirname = os.path.dirname(__file__)
    # rel_path =
    biomartdb = os.path.join(dirname, "data/biomart_GRCh38p13_nov2019.dat")

    with open(biomartdb) as f:
        # get col names
        cols = f.readline().split(',')
        # command line for grep (very fast)
        cmd = 'grep \'' + ensid + '\' ' + biomartdb
        # call shell process
        out, err = call_subprocess(cmd)
        line = out.decode('utf-8')
        # avoid possible errors
        if line != b'':
            # detect wether the input is protein or gene id.
            gene_col = cols.index("Gene stable ID")
            prot_col = cols.index("Protein stable ID\n")
            transcript_col = cols.index("Transcript stable ID")
            # check type of isoform
            isoform_col = cols.index("APPRIS annotation")
            isoform_class = line.split(",")[isoform_col].strip()
            if isoform_filter is None:
                isoform_filter = isoform_class
            # translate ensemblid
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
                elif "ENST" in ensid:
                    protID = line.split(",")[prot_col].strip()
                    geneID = line.split(",")[gene_col].strip()
                    transcriptID = ensid
                    return {'protID': protID, 'geneID': geneID, 'transcriptID': transcriptID}
                # error handling
                else:
                    logger.error(
                        'Input Ensembl ID is neither a protein nor a gene.')
                    raise IOError()
            else:
                logger.error(
                    'Input isoform filter ' + isoform_class + ' does not exist. Please check if you misspelt it.')
                raise IOError()
        # error handling
        logger.error(
            'Could not find Ensembl ID in the data of reference.')
        raise IOError()
