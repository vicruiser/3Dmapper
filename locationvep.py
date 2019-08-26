#!/usr/bin/env python3
#from __future__ import print_function

#VEPfile = "/home/vruizser/PhD/2018-2019/PanCancer_data/vep_output/vep_PanCan_chr_20_1-100000"
import pandas
import os

def FileGenerator(VEPfile):
    
    data = pandas.read_csv(VEPfile, skiprows=42, error_bad_lines=False, sep='\t')

    ids2 = pandas.DataFrame(data.Gene.unique(),  columns=['ids'])

    ids2['VEPfile'] = os.path.basename(VEPfile)

    file_name = "holi.csv"
    
    if os.path.isfile('./holi.csv') :
        ids2.to_csv(file_name,mode='a', sep='\t', header=False, index = False)
    else :
        ids2.to_csv(file_name, sep='\t', header=True, index = False)
        
