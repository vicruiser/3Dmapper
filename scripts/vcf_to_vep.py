#!/usr/bin/python3
import os
import subprocess
import sys

# shell command to execute bcftools
bcftools = "./required_packages/bcftools/bcftools \
+split-vep {} \
-o {} \
-f '%CHROM\_%POS\_%REF\/%ALT \
%CHROM:%POS %Allele %Gene %Feature %Feature_type %Consequence %cDNA_position \
%CDS_position %Protein_position %Amino_acids %Codons %Existing_variation\\n' \
-A tab -d"

add_header = "gawk -i inplace '{if(NR==1){$0=\"#Uploaded_variation Location Allele Gene \
Feature Feature_type Consequence cDNA_position CDS_position Protein_position \
Amino_acids Codons Existing_variation\\n\"$0; print $0}; if(NR!=1){print $0}}' "



# run function
def vcf_to_vep (input_file, out_dir, out_file):
    '''
    Input
    ------

    Output
    ------
    '''
    if os.path.isfile(out_file):
        overwrite = input('File already exists. Overwrite? Y = yes, N = no\n')
        if overwrite.lower() == 'y':
            # export enviroment
            os.environ["BCFTOOLS_PLUGINS"] = "./required_packages/bcftools/plugins"

            # write log file
            log = open(out_dir +'log.txt', 'a')
            log.write('Using bcftools to split the vcf file into vep format...\n')
            log.flush() 

            # command to execute bcftools and add header to file
            cmd = bcftools.format(input_file, out_file)
            cmd2 = add_header + out_file

            # register processes
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log, shell = True)
            p2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell = True)

            # execute processes
            p.communicate()
            p2.communicate()


input_file = '/home/vruizser/PhD/2019-2020/ExAC_benchmark/ExAC.r1.sites.vep.vcf'
input_file = '/home/vruizser/PhD/2019-2020/ExAC_benchmark/head_vcf.txt'
out_dir='./scripts/input_pdbmapper/'
out_file = out_dir + 'converted_vcf.vep2'
vcf_to_vep(input_file, out_dir, out_file)

