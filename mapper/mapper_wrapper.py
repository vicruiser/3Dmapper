# -*- coding: utf-8 -*-
# import necessary modules
import pandas as pd
import numpy as np
import glob
import os
import re
from subprocess import call
import itertools
from .logger import get_logger
from .translate import translate
from .db_parser import parser
from .mapper import mapper
from .decorator import tags
from .run_subprocess import call_subprocess
from .writefile import writefile


# Emojis
DNA = '\U0001F9EC'

# decorator to monitor function
detect_column = "grep -v '##' {} | awk -F ' ' '{{for(i=1;i<=NF;i++) \
{{if ($i ~ /{}/){{print i; exit}}}}}}' "

aa = ['I', 'M','T','N', 'K', 'S', 'R', 'L', 
     'P', 'H', 'Q', 'V','A', 'D',
     'E','G','F', 'Y', 'C', 'W']

# @tags(text_start="Running 3Dmapper...",
#       text_succeed=" Running 3Dmapper...done.",
#       text_fail=" Running 3Dmapper...failed!",
#       emoji=DNA)
def wrapper(id, psdb, vardb, out_dir, pident, evalue, isoform, consequence, loc, index_file, dict_geneprot, varid=None, csv = False, hdf = False):
    
    # logging
    logger = get_logger('wrapper', out_dir)
    # translate ensembl id
    try:
        if id == '-': 
            raise IOError()
        ids = translate(
            id,  out_dir, dict_geneprot, isoform)

        gene_id, prot_id, transcript_id = ids['geneID'], ids['protID'], ids['transcriptID']
        
        if 'APPRIS' in ids.keys():
            APPRIS = ids['APPRIS']
        else:
            APPRIS = list(itertools.repeat(None, len(gene_id))) # change to list with same length as protid
        # run 3Dmapper
        for i in range(0, len(prot_id)):
            try:
                mapper(prot_id[i],
                          gene_id[i],
                          transcript_id[i],
                          psdb,
                          vardb,
                          out_dir,
                          pident,
                          evalue,
                          isoform,
                          APPRIS[i],
                          consequence,
                          loc,
                          varid,
                          csv,
                          hdf)
        # error handling
            except IOError:
                if varid is None:
                    logger.error(
                        ('Warning: {} has no mapping positions.'.format(id)))
                else:
                    logger.error(
                        ('Warning: {} has no mapping positions.'.format(str(varid))))
                continue
        
     # error handling
    except IOError:
        if loc:
            try: 
                if id == '-':
                    transcript_id = '-'
                else: 
                    cmd = ('grep -w \'{}\' {}').format(id, index_file)
                    # call subprocess
                    out, err = call_subprocess(cmd)
                    #print(out, err)
                    if err is None and out != b'':
                        toprocess = out.decode('utf-8')
                        # geneid correspond to the 2nd column in the line
                        cmd2 = detect_column.format(index_file,'Feature')
                        # execute subprocess
                        out2, err2 = call_subprocess(cmd2)
                        # error handling
                        if err2 is None and out2 != b'':
                            col_index_trasncriptid = re.findall('\d+', out2.decode('utf8'))[0]
                            transcript_id = toprocess.split("")[col_index_trasncriptid]
                        #transcript_id = [s for s in toprocess.split(
                        #                    " ") if 'ENST' in s][0]
                annovars_left = parser(transcript_id, vardb)
                
                    #annovars_left = annovars[annovars['Feature']==id]
                            # filter by position type if one or more selected
                if varid is not None:
                    if 'Existing_variation' in annovars_left.columns:
                        annovars_left =  annovars_left[
                            (annovars_left.Uploaded_variation.astype(str) == str(varid)) |
                            ( annovars_left.Existing_variation.astype(str) == str(varid)) ]
                    else:
                        annovars_left =  annovars_left[
                            ( annovars_left.Uploaded_variation.astype(str) == str(varid))]
                    logger.info('position \'' + str(varid) + '\' has been selected.')
                    # if filter returns an empty df, raise error
                    if annovars_left.empty:
                        logger.error(
                            'positions could not be filtered by position id \'' + str(varid) + '\'')
                        raise IOError()
                try: 
                    coding_positions_index = annovars_left.Amino_acids.str.contains('|'.join(aa), regex=True, na = True)
                    noncoding_positions = annovars_left.loc[~coding_positions_index]
                except: 
                    noncoding_positions = False 
                if isoform is None:
                    isoform = ['all']
                if consequence is None:
                    consequence = ['all']  
                # non-protein coding mutations
                if noncoding_positions is not False:
                    noncoding_positions['APPRIS_isoform'] = ''
                    noncoding_positions['Mapping_position'] = 'Noncoding'
                    writefile(transcript_id, out_dir, float(pident), isoform, consequence, noncoding_positions, 'NoncodingPositions', csv, hdf)
                    unmapped_positions = annovars_left.loc[coding_positions_index]
                else: 
                    unmapped_positions = annovars_left
                        
                if unmapped_positions.empty is False:
                    #unmapped_positions = unmapped_positions.iloc[:, 0:16]
                    unmapped_positions.drop_duplicates(inplace=True)
                    unmapped_positions['APPRIS_isoform'] = ''
                    unmapped_positions['Mapping_position'] = 'Unmapped'
                    writefile(transcript_id, out_dir, float(pident), isoform, consequence, unmapped_positions, 'UnmappedPositions', csv, hdf)
            except:
                pass
        logger.error('Warning: {} has no matching ensembl ids.'.format(id))
