#!/usr/bin/python3
# import necessary modules
import os
import os.path
import time
import subprocess
import sys
from scripts.decorator import tags
from timeit import default_timer as timer

# define request function to avoid repeating code


def request(input_file, out_dir, out_file):
    '''
    VCF to VEP format using the plugin "split-vep" from bcftools.

    Parameters
    ----------
    input_file : str        
        Path to infile.
    out_dir : str        
        Path to output.
    out_file : str        
        Name of the output file. 

    Returns
    -------
    converted_vcf.vep 
        Converted file.
    '''
    # export enviroment
    os.environ["BCFTOOLS_PLUGINS"] = "/home/vruizser/PhD/2018-2019/git/PDBmapper/required_packages/bcftools/plugins"

    # shell command to execute bcftools (bash)
    bcftools = "/home/vruizser/PhD/2018-2019/git/PDBmapper/required_packages/bcftools/bcftools \
+split-vep {} \
-o {} \
-f '%CHROM\_%POS\_%REF\/%ALT \
%CHROM:%POS %Allele %Gene %Feature %Feature_type %Consequence %cDNA_position \
%CDS_position %Protein_position %Amino_acids %Codons %Existing_variation\\n' \
-A tab -d"
    # add input variables to command line
    cmd = bcftools.format(input_file, out_file)
    # write log file
    log = open(os.path.join(out_dir, 'log.txt'), 'a')
    log.write('Using bcftools to convert vcf file to vep format...\n')
    log.flush()
    # execute process
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log, shell=True)
    # output from process
    out, err = p.communicate()
    # raise error if process failed
    if err is None:
        pass
    else:
        raise IOError()


# add decorator to main function
@tags(text_start="Converting vcf to vep...This might take up some time...\n",
      text_succeed="Converting vcf to vep...done.\n",
      text_fail="Converting vcf to vep...failed!\n",
      emoji="\U0001F504 ")
def vcf_to_vep(input_file, out_dir, out_file, overwrite):
    '''
    VCF to VEP format using the plugin "split-vep" from bcftools.

    Parameters
    ----------
    input_file : str        
        Path to infile.
    out_dir : str        
        Path to output.
    out_file : str        
        Name of the output file.
    overwrite : str
        Force to overwrite. Default is yes.

    Returns
    -------
    converted_vcf.vep 
        Converted file.
    '''

    # create output dir if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    # execute function
    if os.path.isfile(out_file):
        if overwrite.lower() == 'y':
            request(input_file, out_dir, out_file)
    else:
        request(input_file, out_dir, out_file)
