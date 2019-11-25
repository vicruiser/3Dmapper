#!/bin/sh

#input arguments
VEP_FILE="$1"
OUT_DIR="./input_pdbmapper/splitted_vcf_db/"
mkdir -p $OUT_DIR

echo 'executing'

#detect the column postion containing the protein IDs
COLUMN_INDEX=$(awk -F " " '{for(i=1;i<=NF;i++)\
{if ($i ~ /ENSG/){print i; exit}}}' $VEP_FILE)

echo $COLUMN_INDEX

if [ $COLUMN_INDEX] ; then
        grep -v '##' $VEP_FILE | \
        awk -v ci="$COLUMN_INDEX" \
        -v od="$OUT_DIR" \
        -F ' ' 'NR==1 \
            {h=$0; next} \
            {f=od$ci".vep"} \
            !($ci in p) \
            {p[$ci]; print h > f} \
            {print >> f; close(f)}' 
else
	echo "This file does not contain ENSP ids."
fi


