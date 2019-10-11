#!/bin/bash

#input arguments
VEP_FILE="$1"
OUT_DIR=${2:-./input_pdbmapper/splitted_vcf_db}
mkdir -p $OUTPUT_DIR



VEP_FILE='/home/vruizser/PhD/2018-2019/PanCancer_data/vep_output/vep_PanCan_chr_1_1-100000'
OUT_DIR='./input_pdbmapper/splitted_vcf_db/'
mkdir -p $OUT_DIR


#detect the column postion containing the protein IDs
COLUMN_INDEX=$(awk -F " " '{for(i=1;i<=NF;i++)\
{if ($i ~ /ENSP/){print i; exit}}}'$VEP_FILE)

# split .vep file by protein ID
awk -v ci="$COLUMN_INDEX" \
    -v od="$OUT_DIR" \
    -F ' ' 'NR==1 \
            {h=$0; next} \
            {f=od$ci".vep"} \
            !($ci in p) \
            {p[$ci]; print h > f} \
            {print >> f; close(f)}' \
    $VEP_FILE
