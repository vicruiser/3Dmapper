#!/bin/bash
# if bcftools not installed, install it, 
# activate pluginh to split vep
VCF_FILE="$1"
OUT_FILENAME="$2" #output file name 


#The plugin allows to extract fields from structured annotations such as
# INFO/CSQ created by bcftools/csq or VEP. Assuming the tag added by VEP
# is the INFO/CSQ field, let’s start with printing the list of available
# subfields:

export BCFTOOLS_PLUGINS=./required_packages/bcftools/plugins 
 ./required_packages/bcftools/bcftools +split-vep $VCF_FILE \
                                       -o $OUT_FILENAME.vep \
                                       -f '%CHROM\_%POS\_%REF\/%ALT \
%CHROM:%POS %Allele %Gene %Feature %Feature_type %Consequence %cDNA_position \
%CDS_position %Protein_position %Amino_acids %Codons %Existing_variation \
IMPACT=%IMPACT;DISTANCE=%DISTANCE;STRAND=%STRAND;SIFT=%SIFT\n' -A tab -d

# add header
gawk -i inplace '{if(NR==1){$0="#Uploaded_variation Location Allele Gene \
Feature Feature_type Consequence cDNA_position CDS_position Protein_position \
Amino_acids Codons Existing_variation Extra\n"$0; print $0};
                  if(NR!=1){print $0}}' $OUT_FILENAME.vep

