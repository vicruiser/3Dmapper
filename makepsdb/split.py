# -*- coding: utf-8 -*-
import subprocess
import os
import os.path
import re
import time
from halo import Halo
from .decorator import tags
from .logger import get_logger

detect_column = "awk -F ' ' '{{for(i=1;i<=NF;i++) \
{{if ($i ~ /{}/){{print i; exit}}}}}}' {} "

split_cmd = "awk -v ci=\"{1}\" \
-v od=\"{2}/\" \
-F '\t' 'NR==1 {{h=$0; next}}; \
!seen[$ci]++{{f=od$ci\".{3}\"; print h >> f}}; \
{{f=od$ci\".{3}\"; print >> f; close(f)}}' {0}"

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
    # execute process
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT,
                          shell=True)
    # get output
    out1, err1 = p1.communicate()
    # error handling
    if err1 is None and out1 != b'':
        logger.info('This file contains protein ids')
    else:
        logger.error(err1)
        logger.error(
            'This file could not be splitted. Check the format of your input file.')
        raise IOError()
    # detect if there is output
    col_index = re.findall('\d+', out1.decode('utf8'))[0]
    # stop if no ENSP id detected
    if col_index != '':
        # write log file
        logger.info('Splitting interfaces file...')
        # Second command
        if parallel is True:
            cmd2 = split_cmd_parallel.format(input_file, col_index,
                                             out_dir, out_extension)
        else:
            cmd2 = split_cmd.format(input_file, col_index,
                                    out_dir, out_extension)
        # register process
        p2 = subprocess.Popen(cmd2,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,
                              shell=True)
        # error handling
        out2, err2 = p2.communicate()
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
