#!/usr/bin/python3
import os
import subprocess
import sys


def request(input_file, out_dir, out_file):
    # export enviroment
    os.environ["BCFTOOLS_PLUGINS"] = "/home/vruizser/PhD/2018-2019/git/PDBmapper/required_packages/bcftools/plugins"
    
    # shell command to execute bcftools
    bcftools = "/home/vruizser/PhD/2018-2019/git/PDBmapper/required_packages/bcftools/bcftools \
+split-vep {} \
-o {} \
-f '%CHROM\_%POS\_%REF\/%ALT \
%CHROM:%POS %Allele %Gene %Feature %Feature_type %Consequence %cDNA_position \
%CDS_position %Protein_position %Amino_acids %Codons %Existing_variation\\n' \
-A tab -d"

    # write log file
    log = open(out_dir +'log.txt', 'a')
    log.write('Using bcftools to split the vcf file into vep format...\n')
    log.flush() 

    # register processes
    cmd = bcftools.format(input_file, out_file)
    subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log, shell = True)


# run function
def vcf_to_vep (input_file, out_dir, out_file, overwrite):
    '''
    Input
    ------

    Output
    ------
    '''
    #create dir if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)

    if os.path.isfile(out_file):

        if overwrite.lower() == 'y':
            request(input_file, out_dir, out_file)

        else:
            print('Using already created converted_vcf.vep file.')

    else:
        request(input_file, out_dir, out_file)


