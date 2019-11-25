#!/usr/bin/python3
import subprocess
import os
import os.path
import re
import time
from halo import Halo
from decorator import tags

detect_column = "awk -F ' ' '{{for(i=1;i<=NF;i++) \
{{if ($i ~ /{}/){{print i; exit}}}}}}' {} "

split_cmd = "grep -v '##' {} | \
sed -e '1s/^#//' | \
awk -v ci=\"{}\" \
-v od=\"{}/\" \
-F ' ' 'NR==1 \
{{h=$0; next}} \
{{f=od$ci\".{}\"}} !($ci in p) \
{{p[$ci]; print h > f}} \
{{print >> f; close(f)}}'"


def request(prefix, input_file, out_dir, out_extension):
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
    log1 = open(os.path.join(out_dir, 'log_find_column.txt'), 'a')
    log1.write('Find column...\n')
    log1.flush()
    # First command
    cmd1 = detect_column.format(prefix, input_file)
    # execute process
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE,
                          shell=True)
    # get output
    out1, err1 = p1.communicate()
    # detect if there is output
    col_index = re.findall('\d+', out1.decode('utf8'))[0]
    # stop if no ENSP id detected
    if col_index != '':
        # write log file
        log2 = open(os.path.join(out_dir, 'log_split_files.txt'), 'a')
        log2.write('Splitting file...\n')
        log2.flush()
        # Second command
        cmd2 = split_cmd.format(input_file, col_index, out_dir, out_extension)
        # register process
        p2 = subprocess.Popen(cmd2,
                              stdout=subprocess.PIPE,
                              stderr=log2,
                              shell=True)
        # error handling
        out2, err2 = p2.communicate()
        if err2 is not None:
            raise IOError()
    else:
        raise IOError()


# add decorator to main function
@tags(text_start="Split file by selected ensembl id...This might take up some time...",
      text_succeed="Split file by selected ensembl id...done.",
      text_fail="Split file by selected ensembl id...failed!",
      emoji="\U00002702")
def split(prefix, input_file, out_dir, out_extension, overwrite):
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
            request(prefix, input_file, out_dir, out_extension)
    else:
        request(prefix, input_file, out_dir, out_extension)
