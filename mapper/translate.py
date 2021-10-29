# -*- coding: utf-8 -*-

# import necessary modules
import sys
import os
import re
import pandas as pd
import numpy as np
from .run_subprocess import call_subprocess
from .logger import get_logger


def translate(id, log_dir,dict_geneprot=None, isoform_filter=None):
    '''
    gene-transcript-protein id translator.

    Parameters
    ----------
    id : str
        Input ID corresponding to a gene, transcript or protein.

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
    #if dict_geneprot is None : 
    #dict_geneprot = os.path.join(dirname, "data/biomart_GRCh38p13_feb2021.dat")

    with open(dict_geneprot) as f:
        # get col names
        cols = f.readline().rstrip().split(',')
        # command line for grep (very fast)
        cmd = 'grep -w \'' + id + '\' ' + dict_geneprot
        # call shell process
        out, err = call_subprocess(cmd)
        line = out.decode('utf-8')
        # avoid possible errors
        if out != b'':
            list_ids = line.split('\n')[:-1]
            df = pd.DataFrame(np.array(list_ids), columns=list("l"))
            #df[['isoform', 'geneID', 'transcriptID', 'protID', 'UniprotID']
            #   ] = df['l'].str.split(',', expand=True)
            df[cols
               ] = df['l'].str.split(',', expand=True)
            # filter by principal isoform if any filter
            protID = df['protID'].tolist()
            geneID = df['geneID'].tolist()
            transcriptID = df['transcriptID'].tolist()
            results = {'protID': protID, 'geneID': geneID, 'transcriptID': transcriptID}
            if isoform_filter is not None:
                df = df[df['isoform'].isin(isoform_filter)]
                if df.empty:
                    logger.error(
                        'Input isoform filter ' + isoform_filter + ' does not exist. Please check if you misspelt it.')
                    raise IOError()
                else : 
                    APPRIS = df['isoform'].tolist()
                    results['APPRIS'] = APPRIS
                  
            return results

        else:
            logger.error(
                'Input Ensembl ID is neither a protein nor a gene.')
            raise IOError
            
