#!/bin/bash

VCF_FILE="$1"
OUT_FILENAME="$2"

# if detect INFO column or 6 columns?
echo 'Your input is a BCF file. Transforming to VEP format...'
sh vcf_to_vep.sh $VCF_FILE $OUT_FILENAME
echo 'File successfuly converted.'

# split always
echo 'Splitting VEP file into  individual protein IDs.'
sh split_vep_by_protid.sh $VEP_file
echo 'Done.'


# now we can run PDBmapper
#python3 main.py --protid ENSP00000482258 -vcf 