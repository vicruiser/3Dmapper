#!/usr/bin/python3
import subprocess
import os

detect_column = "awk -F ' ' '{{for(i=1;i<=NF;i++) \
{{if ($i ~ /ENSP/){{print i; exit}}}}}}' {} "

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
                        
    # Second command
    cmd2 = split.format(input_file, out.decode('utf8'), out_dir) 
    #registe process
    p2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell = True)
    # Execute processes
    p2.communicate()

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

                    


#input_file = ['/home/vruizser/PhD/2018-2019/PanCancer_data/vep_output/vep_PanCan_chr_1_1-100000',
#              '/home/vruizser/PhD/2018-2019/PanCancer_data/vep_output/vep_PanCan_chr_1_100001-200000']

#out_dir='./scripts/input_pdbmapper/splitted_vcf_db/'
#print("heloo")

#for f in input_file: 
#    split_vep(f, out_dir) # overwrite !!!!
