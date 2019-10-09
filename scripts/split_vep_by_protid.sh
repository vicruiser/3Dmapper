#!/bin/bash

#input arguments
VEP_FILE="$1"
mkdir -p ./input/splitted_vep_db
OUTPUT_DIR='./input/splitted_vep_db/'

#detect the column postion containing the protein IDs
COLUMN_INDEX=$(awk -F " " '{for(i=1;i<=NF;i++)\
                            {if ($i ~ /ENSP/)\
                                {print i; exit}\
                                }\
                            }'\
                    $VEP_FILE)

# split .vep file by protein ID
awk -v ci="$COLUMN_INDEX" \
    -v od="$OUTPUT_DIR" \
    -F ' ' 'NR==1 \
            {h=$0; next} \
            {f=od$ci".vep"} \
            !($ci in p) \
            {p[$ci]; print h > f} \
            {print >> f; close(f)}' \
    $VEP_FILE
