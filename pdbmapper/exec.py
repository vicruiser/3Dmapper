# -*- coding: utf-8 -*-
# import necessary modules
import os


def exec(ensemblid, filter_iso, out, *args, **kwargs, args.intdb,
         args.vardb,
         args.out,
         args.pident,
         args.filter_var, log_dir):

    try:
        ensemblIDs = translate_ensembl(
            ensemblid, args.filter_iso, args.out)
        geneid = ensemblIDs['geneID']
         protid = ensemblIDs['protID']
          transcriptID = ensemblIDs['transcriptID']
           # run PDBmapper
           try:
                PDBmapper(protid,
                          geneid,
                          transcriptID,
                          args.intdb,
                          args.vardb,
                          args.out,
                          args.pident,
                          args.filter_var)
        # error handling
            except IOError:
                logger.error('Warning: ' + ensemblid +
                             ' has no mapping variants.')
