# -*- coding: utf-8 -*-

# import necessary modules
import sys
import os
import re
import pandas as pd
import numpy as np
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
        #cols = f.readline().split(',')
        # command line for grep (very fast)
        cmd = 'grep \'' + ensid + '\' ' + biomartdb
        # call shell process
        out, err = call_subprocess(cmd)
        line = out.decode('utf-8')
        # avoid possible errors
        if line != b'':
            list_ids = line.split('\n')
            df = pd.DataFrame(np.array(list_ids), columns=list("l"))
            df[['isoform', 'geneID', 'transcriptID', 'protID']
               ] = df['l'].str.split(',', expand=True)

            # filter by principal isoform if any filter
            if isoform_filter is not None:
                df = df[df['isoform'].isin(isoform_filter)]
                if df.empty:
                    logger.error(
                        'Input isoform filter ' + isoform_class + ' does not exist. Please check if you misspelt it.')
                    raise IOError()

            protID = df['protID'].tolist()
            geneID = df['geneID'].tolist()
            transcriptID = df['transcriptID'].tolist()
            return {'protID': protID, 'geneID': geneID, 'transcriptID': transcriptID}

        else:
            logger.error(
                'Input Ensembl ID is neither a protein nor a gene.')
            raise IOError()
            # else:
            #     logger.error(
            #         'Input isoform filter ' + isoform_class + ' does not exist. Please check if you misspelt it.')
            #     raise IOError()
        # error handling
        # logger.error(
        #    'Could not find Ensembl ID in the data of reference.')
        #raise IOError()
