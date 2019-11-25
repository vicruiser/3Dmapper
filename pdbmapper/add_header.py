# -*- coding: utf-8 -*-
import os
import subprocess
import sys
from decorator import tags


def add_header(vep_file):
    '''
    Add header to vep file generated with vcf_to_vep.py

    Parameters
    ----------
    vep_file : str
        input file
    Returns
    -------
    file with header

    '''
    #
    add_header = "gawk -i inplace '{{if(NR==1){{$0=\"Uploaded_variation Location Allele Gene \
Feature Feature_type Consequence cDNA_position CDS_position Protein_position \
Amino_acids Codons Existing_variation\\n\"$0; print $0}}; if(NR!=1){{print $0}}}}' {}"

    cmd2 = add_header.format(vep_file)
    p2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)
    p2.wait()
