#!/usr/bin/python3
import subprocess
import os
import re
import threading
import subprocess
#import itertools, sys
import time
from halo import Halo

spinner = Halo(text='Loading', spinner='dots')



detect_column = "awk -F ' ' '{{for(i=1;i<=NF;i++) \
{{if ($i ~ /{}/){{print i; exit}}}}}}' {} "

split_cmd = "grep -v '##' {} | \
awk -v ci=\"{}\" \
-v od=\"{}\" \
-F ' ' 'NR==1 \
{{h=$0; next}} \
{{f=od$ci\".{}\"}} !($ci in p) \
{{p[$ci]; print h > f}} \
{{print >> f; close(f)}}'"


def request(prefix, input_file, out_dir, out_extension):
    # First command
    cmd = detect_column.format(prefix, input_file)
    # register process
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
    # execute process
    out, err = p1.communicate()
    n = re.findall('\d+', out.decode('utf8'))[0]
    

    # stop if no ENSP id detected
    if n != '' :
        # Second command
        cmd2 = split_cmd.format(input_file, n, out_dir, out_extension)
        # write log file
        log = open(out_dir +'log.txt', 'a')
        log.write('Using bcftools to split the vcf file into vep format...\n')
        log.flush() 

        #register process
        p = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr = log, shell = True)

        while p.poll() is None:
            spinner.start( "Splitting file")

        spinner.stop()
        spinner.succeed("Now you have your database")

        out1, err1 = p.communicate()
        print(out1.decode('utf-8'))
        if(err1 is not None):
	        print(err1.decode('utf-8'))    
    else: 
        spinner.fail("Failed")
        raise IOError()

# funciton
def split (input_file, out_dir, overwrite, prefix, out_extension):
    '''
    Input
    ------

    Output
    ------
    '''
    #create dir if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    # check if this process has been already executed. 
    if any(f.endswith("." + out_extension) for f in os.listdir(out_dir)):
        if overwrite.lower() == 'y':    
            print("hola")
            try:
                request(prefix, input_file, out_dir, out_extension)
            except IOError:
                print("que ha pachao")
                spinner.stop( "putamerda")

    else: 
        print("hola")
        request(prefix, input_file, out_dir, out_extension)
