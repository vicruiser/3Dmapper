# -*- coding: utf-8 -*-
# import necessary modules
import os
import os.path
import time
import subprocess
import sys
from .decorator import tags
from timeit import default_timer as timer
from .logger import get_logger


# define request function to avoid repeating code
def request(input_file, out_dir, out_file, log_dir, parallel=False):
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
    os.environ["BCFTOOLS_PLUGINS"] = "%s/bin/" % os.environ.get(
        'VIRTUAL_ENV', '/usr/local/')

    # shell command to execute bcftools (bash)
    bcftools = "bcftools \
+split-vep {} \
-o {} \
-f '%CHROM\_%POS\_%REF\/%ALT \
%CHROM:%POS %Allele %Gene %Feature %Feature_type %Consequence %cDNA_position \
%CDS_position %Protein_position %Amino_acids %Codons %Existing_variation\\n' \
-A tab -d"

    bcftools_parallel = "parallel --pipe -q bcftools \
+split-vep {} \
-o {} \
-f '%CHROM\_%POS\_%REF\/%ALT \
%CHROM:%POS %Allele %Gene %Feature %Feature_type %Consequence %cDNA_position \
%CDS_position %Protein_position %Amino_acids %Codons %Existing_variation\\n' \
-A tab -d"
    # add input variables to command line
    if parallel is True:
        cmd = bcftools_parallel.format(input_file, out_file)
    else:
        cmd = bcftools.format(input_file, out_file)
    # log file
    logger = get_logger('vcf2vep', log_dir)
    logger.info('Using bcftools to convert vcf file to vep format.')

    # execute subprocess
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, shell=True)
    # output from process
    out, err = p.communicate()

    # raise error if process failed
    if err is None:
        logger.info('File in vcf format converted successfully to vep format.')
    else:
        logger.error(err)
        logger.error(
            'Something went wrong converting the vcf file into vep format.')
        raise IOError()


# add decorator to main function
@tags(text_start="Converting vcf to vep...This might take up some time...\n",
      text_succeed="Converting vcf to vep...done.\n",
      text_fail="Converting vcf to vep...failed!\n",
      emoji="\U0001F504 ")
def vcf2vep(input_file, out_dir, out_file, overwrite, log_dir, parallel=False):
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
        if overwrite is True:
            request(input_file, out_dir, out_file, overwrite, log_dir, parallel)
    else:
        request(input_file, out_dir, out_file,  overwrite, log_dir, parallel)
