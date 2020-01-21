# -*- coding: utf-8 -*-
# import necessary modules
import os
from .logger import get_logger
from .translate_ensembl import translate_ensembl
from .PDBmapper import PDBmapper
from .decorator import tags

# Emojis
DNA = '\U0001F9EC'

# decorator to monitor function


# @tags(text_start="Running PDBmapper...",
#       text_succeed=" Running PDBmapper...done.",
#       text_fail=" Running PDBmapper...failed!",
#       emoji=DNA)
def wrapper(ensemblid, intdb, vardb, out, pident, feature_type, varid=None):

    # logging
    logger = get_logger('wrapper', out)
    # translate ensembl id
    try:
        ensemblIDs = translate_ensembl(
            ensemblid,  out, feature_type)
        geneid = ensemblIDs['geneID']
        protid = ensemblIDs['protID']
        transcriptID = ensemblIDs['transcriptID']
     # error handling
    except:
        logger.error(
            ('Warning: {} has no matching ensembl ids.'.format(ensemblid)))
    # run PDBmapper
    try:
        PDBmapper(protid,
                  geneid,
                  transcriptID,
                  intdb,
                  vardb,
                  out,
                  pident,
                  feature_type,
                  varid)
    # error handling
    except:
        if varid is None:
            logger.error(
                ('Warning: {} has no mapping variants.'.format(ensemblid)))
        else:
            logger.error(
                ('Warning: {} has no mapping variants.'.format(str(varid))))
