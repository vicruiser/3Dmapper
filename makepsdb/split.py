# -*- coding: utf-8 -*-
import subprocess
import os
import os.path
import re
import time
from halo import Halo
from .decorator import tags
from .logger import get_logger
from .run_subprocess import call_subprocess

import pandas as pd
import glob
import os


#sort_cmd="awk 'NR==1; NR>1{print $0 | \"sort -n\"}' headd.txt > head_sorted.txt"

detect_column = "awk -F '\t' '{{for(i=1;i<=NF;i++) \
{{if ($i ~ /{}/){{print i; exit}}}}}}' {} "

split_cmd = "awk -v ci=\"{1}\" \
-v od=\"{2}/\" \
-F '\t' 'NR==1 {{h=$0; next}}; \
!seen[$ci]++{{f=od$ci\".{3}\"; print h > f}}; \
{{f=od$ci\".{3}\"; print >> f}}' {0}"
#; close(f)
#split_cmd= "awk -v c=\"{1}\" -v od=\"{2}/\" \
#-F'\t' \
#'NR==1{h=$0; next} \
#!seen[$c]++{print h > od$c\".txt\"} \
#{print >> od$c\".txt\"}' {0}"

split_cmd_parallel = "gawk -v ci=\"{1}\" \
-v od=\"{2}/\" \
-F '\t' 'NR==1 {{h=$0; next}}; \
!seen[$ci]++{{f=od$ci\".{3}\"; print h >> f}}; \
{{f=od$ci\".{3}\"; print >> f; close(f)}}' {0}"

def request(prefix, input_file, out_dir, out_extension, log_dir, parallel=False):
    '''
    VCF to VEP format using the plugin "split-vep" from bcftools.

    Parameters
    ----------
    prefix : str
        Ensembl ID prefix. Either ESNG or ENSP.
    input_file : str
        Path to infile.
    out_dir : str
        Path to output.
    out_extension : str
        Output filename extension

    Returns
    -------
    ./dir
        Directory containing splitted files.
    '''
    # log file
    logger = get_logger('split', log_dir)
    logger.info('Splitting input file.')
    # First command
    cmd1 = detect_column.format(prefix, input_file)
    
    # command
   # cmd1 = detect_column.format(input_file, prefix)
    # execute subprocess
    out1, err1 = call_subprocess(cmd1)
    # execute process
    #p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE,
    #                      stderr=subprocess.STDOUT,
    #                      shell=True)
    # get output
    #out1, err1 = p1.communicate()
    # error handling
    if err1 is None and out1 != b'':
        col_index = re.findall('\d+', out1.decode('utf8'))[0]
        logger.info('This file contains protein ids')
        # detect if there is output
    else:
        logger.error(err1)
        logger.error(
            'This file could not be splitted. Check the format of your input file.')
        raise IOError()
 
    # stop if no ENSP id detected
    if col_index != '':
        # write log file
        logger.info('Splitting interfaces file...')

        # maybe play around to get better performance 
        #chunksize = 10000000

        #files = glob.glob('./file_*.csv')
        

        #for chunk in pd.read_csv(input_file, chunksize=chunksize,  sep="\t" ):
            #print(chunk)
        #    u_inst = chunk['Protein_accession'].unique()
        #    for inst in u_inst:
                # filter instrument data
        #        inst_df = chunk[chunk.Protein_accession == inst]
                # filter columns
                #inst_df = inst_df[['time', 'code', 'val']]
                # append to instrument file
                # only write header if not exist yet
        #        inst_file = os.path.join(out_dir, inst +'.txt')
        #        file_exist = os.path.isfile(inst_file)
        #        inst_df.to_csv(inst_file, mode='a', header=not file_exist, sep = '\t', index = False)

        # Second command
        # f = open(input_file, 'r')       
        # header = f.readline()
        # #pid1 = ''
        # for line in f: 
        #     #print(line)
        #     pid = line.split('\t')[int(col_index)-1]
        #     if os.path.isfile(os.path.join(out_dir, pid +'.txt')) is False: 
        #         with open(os.path.join(out_dir, pid +'.txt'), 'a') as outfile: 
        #             outfile.write(header)
        #     #if pid1 != pid2: 
        #     else:
        #         with open(os.path.join(out_dir, pid +'.txt'), 'a') as outfile: 
        #             outfile.write(line)
            #else : 
            #    with open(os.path.join(out_dir, pid2+'.txt'), 'a') as outfile: 
            #        outfile.write(line)
            #pid1 = pid2

            
            
        
        
        if parallel is True:
            cmd2 = split_cmd_parallel.format(input_file, col_index,
                                             out_dir, out_extension)
        else:
            cmd2 = split_cmd.format(input_file, col_index,
                                    out_dir, out_extension)
        # register process
        out2, err2 = call_subprocess(cmd2)
   
        if err2 is None:
            logger.info('This file was splitted successfully.')
        else:
            logger.error(
                'This file could not be splitted. Check the format of your input file.')
            raise IOError()
    else:
        logger.error('The input file has zero protein entries.')
        raise IOError()


# add decorator to main function
@tags(text_start="Creating protein structures DB...This might take up some time...",
      text_succeed="Creating protein structures DB...done.",
      text_fail="Creating protein structures DB...failed!. Check the format of your input file.",
      emoji="\U00002702")
def split(prefix, input_file, out_dir, out_extension, overwrite, log_dir, parallel=False):
    '''
    VCF to VEP format using the plugin "split-vep" from bcftools.

    Parameters
    ----------
    prefix : str
        Ensembl ID prefix. Either ESNG or ENSP. 
    input_file : str        
        Path to infile.
    out_dir : str        
        Path to output.
    out_extension : str        
        Output filename extension 
    overwrite : str
        Force to overwrite. Default is yes.

    Returns
    -------
    ./dir 
        Directory containing splitted files. 
    '''
    # create dir if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    # execute request function
    if any(f.endswith("." + out_extension) for f in os.listdir(out_dir)):

        if overwrite is True:
            request(prefix, input_file, out_dir,
                    out_extension, log_dir, parallel)
    else:
        request(prefix, input_file, out_dir, out_extension, log_dir, parallel)
