#!/usr/bin/python3
import subprocess
import os

detect_column = "awk -F ' ' '{{for(i=1;i<=NF;i++) \
{{if ($i ~ /ENSG/){{print i; exit}}}}}}' {} "

split = "grep -v '##' {} | \
awk -v ci=\"{}\" \
-v od=\"{}\" \
-F ' ' 'NR==1 \
{{h=$0; next}} \
{{f=od$ci\".vep\"}} !($ci in p) \
{{p[$ci]; print h > f}} \
{{print >> f; close(f)}}'"


def request(input_file, out_dir):
    # First command
    cmd = detect_column.format(input_file)
    # register process
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
    # execute process
    out, err = p1.communicate()
    # stop if no ENSP id detected
    if out.decode('utf8') != '' :
        # Second command
        cmd2 = split.format(input_file, out.decode('utf8'), out_dir) 
        #registe process
        p2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell = True)
        p2.wait()
    else: 
        raise IOError()

# funciton
def split_vep (input_file, out_dir, overwrite):
    '''
    Input
    ------

    Output
    ------
    '''
    #create dir if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)

    # check if this process has been already executed. 
    if any(f.endswith(".vep") for f in os.listdir(out_dir)):
        if overwrite.lower() == 'y':    
            request(input_file, out_dir)

    else: 
        request(input_file, out_dir)
