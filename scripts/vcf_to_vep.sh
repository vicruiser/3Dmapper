#!/bin/bash
# if bcftools not installed, install it, 
# activate pluginh to split vep
VCFfile="$1"



#The plugin allows to extract fields from structured annotations such as INFO/CSQ created by bcftools/csq or VEP.
#Assuming the tag added by VEP is the INFO/CSQ field, let’s start with printing the list of available subfields:

if which programname >/dev/null; then
    bcftools +split-vep $VCFfile -f '%CHROM:%POS %CSQ\n' -d -A tab
    echo exists
else
    echo does not exist
fi

