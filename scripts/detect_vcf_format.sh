#!/bin/bash
VCF_FILE="$1"

if grep -q CSQ $VCF_FILE; then
    return 'VCF'
else
    return 'VEP'
fi