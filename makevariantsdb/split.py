# -*- coding: utf-8 -*-
import subprocess
import os
import os.path
import re
import time
import glob
from halo import Halo
from .decorator import tags
from .logger import get_logger
from .parse_argv import parse_commandline
#from .create_var_index import index


# for all the columns, find the one that matches with the pattern ENSG
detect_column = "awk -F ' ' '{{for(i=1;i<=NF;i++) \
{{if ($i ~ /{}/){{print i; exit}}}}}}' {} "

split_cmd = "grep -v '##' {} | \
sed -e '1s/^#//' | \
awk -v ci=\"{}\" \
-v od=\"{}/\" \
-F ' ' 'NR==1 {{h=$0; next}} \
{{f=od$ci\".{}\"}} \
!($ci in p) {{p[$ci]}} \
system(\" stat \" f \" > /dev/null 2> /dev/null\") != 0 {{print h > f }} \
{{print >> f; close(f)}}'"

index_file = "awk -F ' ' 'NR>2{{print ${}, ${}, ${}}}' {} | uniq >> {}  "

# - grep -v '##': remove lines starting with ##
# - sed -e '1s/^#// : remove # from header line
# - awk
#       -v ci=\"{}\": set variable ci
#       -v od=\"{}/\" : set variable od
#       -F ' ' : file separated by ' '
#       - 'NR==1 {{h=$0; next}} : first line is header (h)
#       - {{f=od$ci\".{}\"}} : set variable f as output filename
#       - !($ci in p) {{p[$ci]}} : if  number of column index in file, select it
#       - system(\"stat \" f \" >/dev/null 2>/dev/null\") != 0 {{print h > f }} : add header to outfile only if out file doesn't exist
#       - {{print >> f; close(f)}}'" print line to file


def request(prefix, input_file, out_dir, out_extension, log_dir):
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
    if err1 is None:
        logger.info('This file contains gene ids')
    else:
        logger.error(err1)
        logger.error('This file could not be splitted')
        raise IOError()
    
    # First command
    cmd2 = detect_column.format('Uploaded_variation', input_file)
    # execute process
    p2 = subprocess.Popen(cmd2,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT,
                          shell=True)
    # get output
    out2, err2 = p2.communicate()
    # error handling
    if err2 is None:
        logger.info('Indexing variants file.')
    else:
        logger.error(err2)
        logger.error('This file cannot be indexed. \
            Does not contain a column with variants ids.')
        raise IOError()
    
    # First command
    cmd3 = detect_column.format('Existing_variation', input_file)
    # execute process
    p3 = subprocess.Popen(cmd3,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT,
                          shell=True)
    # get output
    out3, err3 = p3.communicate()
    # error handling
    if err3 is None:
        logger.info('Indexing variants file.')
    else:
        logger.error(err3)
        logger.error('This file cannot be indexed. \
            Does not contain a column with variants ids.')
        raise IOError()

    # detect if there is output
    col_index_geneid = re.findall('\d+', out1.decode('utf8'))[0]
    col_index_varid = re.findall('\d+', out2.decode('utf8'))[0]
    col_index_namevarid = re.findall('\d+', out3.decode('utf8'))[0]
    # stop if no ENSG id detected
    if col_index_geneid != '':
        # Second command
        cmd4 = split_cmd.format(input_file, col_index_geneid, out_dir, out_extension)
        # register process
        p4 = subprocess.Popen(cmd4,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,
                              shell=True)
        # error handling
        out4, err4 = p4.communicate()

        if err4 is None:
            logger.info('This file was splitted successfully')
        else:
            logger. error(err4)
            logger.error('This file could not be splitted')
            raise IOError()
        
        cmd5=index_file.format(col_index_varid,
                               col_index_geneid,
                               col_index_namevarid,
                               input_file,
                               os.path.join(out_dir,'variants.index'))
        # register process
        p5 = subprocess.Popen(cmd5,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,
                              shell=True)
        # error handling
        out5, err5 = p5.communicate()

        if err5 is None:
            logger.info('This file was indexed.')
        else:
            logger.error(err5)
            logger.error('This file could not be indexed')
            raise IOError()
    else:
        logger.error('The input file has zero gene entries.')
        raise IOError()


# add decorator to main function
@tags(text_start="Split file by selected ensembl id...This might take up some time...",
      text_succeed="Split file by selected ensembl id...done.",
      text_fail="Split file by selected ensembl id...failed!",
      emoji="\U00002702")
def split(prefix, input_file, out_dir, out_extension, overwrite, log_dir):
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

        if overwrite.lower() == 'y':
            request(prefix, input_file, out_dir, out_extension, log_dir)
            
    else:
        request(prefix, input_file, out_dir, out_extension, log_dir)
        
