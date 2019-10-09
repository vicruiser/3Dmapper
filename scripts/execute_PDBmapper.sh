#!/bin/bash

VCF_FILE="$1"

# if detect INFO column or 6 columns?
echo 'Your input is a BCF file. Transforming to VEP format...'
sh vcf_to_vep.sh $VCF_FILE
echo 'File successfuly converted.'

# split always
echo 'Splitting VCF file by individual protein IDs.'
sh split_vep_by_protid.sh $VEP_file
echo 'Done.'


# now we can run PDBmapper
python3 PDBmapper --protid ENSP00000482258 -vcf 