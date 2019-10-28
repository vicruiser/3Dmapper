#!/usr/bin/python3
import os
import subprocess
import sys
from scripts.decorator import tags

def add_header(out_file):
    # add header to resulting vep file
    add_header = "gawk -i inplace '{{if(NR==1){{$0=\"#Uploaded_variation Location Allele Gene \
Feature Feature_type Consequence cDNA_position CDS_position Protein_position \
Amino_acids Codons Existing_variation\\n\"$0; print $0}}; if(NR!=1){{print $0}}}}' {}"
    
    cmd2 = add_header.format(out_file)
    p2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True)
    p2.wait()
