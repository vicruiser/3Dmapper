#!/usr/bin/env python3
# coding: utf-8

# Import necesary modules
import numpy as np
import sys
import os
import gzip
import re
import pandas as pd
from timeit import default_timer as timer
import glob
from subprocess import call
import vcfpy

# import functions from scripts
from scripts.parse_argv import parse_commandline
from scripts.VEPcrossref import VEPfileCrossrefGenerator as cr
import scripts.mapping_tools as mt
import scripts.PDBmapper
import scripts.vep_parser
from scripts.detect_vcf_format import detect_format
from scripts.vcf_to_vep import vcf_to_vep
from scripts.split_vep import split_vep

# define main function to execute the previous defined functions together
def main():
    
    # parse command line options
    args = parse_commandline()
    description = '''

 ------------------------------------------------------------------------------------------------------------------ 

        $$$$$$$\  $$$$$$$\  $$$$$$$\                                                                
        $$  __$$\ $$  __$$\ $$  __$$\                                                                 
        $$ |  $$ |$$ |  $$ |$$ |  $$ |$$$$$$\$$$$\   $$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\  
        $$$$$$$  |$$ |  $$ |$$$$$$$\ |$$  _$$  _$$\  \____$$\ $$  __$$\ $$  __$$\ $$  __$$\ $$  __$$\ 
        $$  ____/ $$ |  $$ |$$  __$$\ $$ / $$ / $$ | $$$$$$$ |$$ /  $$ |$$ /  $$ |$$$$$$$$ |$$ |  \__|
        $$ |      $$ |  $$ |$$ |  $$ |$$ | $$ | $$ |$$  __$$ |$$ |  $$ |$$ |  $$ |$$   ____|$$ |      
        $$ |      $$$$$$$  |$$$$$$$  |$$ | $$ | $$ |\$$$$$$$ |$$$$$$$  |$$$$$$$  |\$$$$$$$\ $$ |      
        \__|      \_______/ \_______/ \__| \__| \__| \_______|$$  ____/ $$  ____/  \_______|\__|      
                                                              $$ |      $$ |                          
                                                              $$ |      $$ |                          
                                                              \__|      \__|   
        
        
------------------------  Map annotated genomic variants to protein interfaces data in 3D. ------------------------

'''
    print(description)  

    # create output directory if it doesn't exist
    if not os.path.exists(args.out):
        os.mkdir(args.out)
        print("Directory", args.out,  "created.")
    else:    
        print(args.out,
              "is an existing directory. Results will be written in there.")
    
    # get geneID
    #biomartdb = open('./default_input_data/gene_transcript_protein_ens_ids.txt', 'r')
    #print('Biomart file read...')
    #ensemblIDs = mt.ensemblID_translator(biomartdb, args.protid)
    #geneID = ensemblIDs['geneID']

     # variants input:
    if args.vep:
        print('''VEP option will be available soon. Using VarMap db instead.
        Otherwise, please provide your own vcf file with the -vcf option.''')
        try:
            # ./vep -i input.vcf -o out.txt -offline
            print('VarMap db is used.')
            ClinVarDB = './default_input_data/splitted_ClinVar'
            VCF_subset = mt.parser(ClinVarDB, geneID, "\t", 'varmap')
        except:
            print('ERROR: cannot open or read input VarMap db file.')
            exit(-1)

    elif args.vcf:
        try:
            print('Detecting inputs variants file format...')
            # for loop in case we have multiple inputs to read
            if len(args.vcf) == 1:
                args.vcf= [args.vcf]
            for f in args.vcf: 
                print(f)
                print(args.vcf)
                # detect the format of the vcf file(s)
                # possible formats: .vcf or .vep
                input_format = detect_format(f)
                if input_format == "vcf":
                    # 1) put vcf into vep format
                    out_dir = args.out + '/pdbmapper/input/' #created by default
                    out_file = out_dir + 'converted_vcf.vep'
                    vcf_to_vep(f, out_dir, out_file, args.force)
                    # 2) split vep file by protein id to speed up the
                    # mapping proccess
                    vcf_db_dir = out_dir + 'vcf_db/' #created by default
                    split_vep(f, vcf_db_dir, args.force)
                elif input_format == "vep":
                    # split vep file by protein id to speed up the
                    # mapping proccess
                    vcf_db_dir = './pdbmapper/input/vcf_db/'
                    split_vep(f, vcf_db_dir, args.force) 
                else: 
                    print(input_format)
        except IOError:
            print('ERROR: cannot open or read input vcf file.')
            exit(-1)
    elif args.varmap:
        try:
            print('VarMap db is used.')
            ClinVarDB = './dbs/splitted_ClinVar/'
            VCF_subset = mt.parser(ClinVarDB, geneID, "\t", 'varmap')
        except IOError:
            print('ERROR: cannot open or read input VarMap db file.')
            exit(-1)

#     # set interfacesDB_susbet:
#     if args.intdb:
#         try:
#             interfacesDB = open(args.intdb, 'r')
#         except IOError:
#             print('ERROR: cannot open or read input interfaces db file.')
#             exit(-1)
#     else:
#         print('Default interfaces DB is used.')
#         interfaces_dir = './dbs/splitted_interfaces_db'
#         interfacesDB_subset = interfaceParse(interfaces_dir, args.protid)


#     # create chimera scripts:
    if args.chimera is not None:
        #chimera()
        pass

#     # run PDBmapper
#     if args.protid:
#         # input list of proteins
#         if isinstance(args.protid, list):
#             # iterate over list
#             for one_protid in args.protid:
#                 try:
#                     PDBmapper(one_protid, interfacesDB_subset, VCF_subset,
#                             args.out)
#                 except IOError:
#                     print('ERROR: Ensembl protein id provided is no supported.')
#                     exit(-1)
#         # single protein id as input
#         else:
#             try:
#                 PDBmapper(args.protid, interfacesDB_subset, VCF_subset,
#                           args.out)
#             except IOError:
#                 print('ERROR: Ensembl protein id provided is no supported.')
#                 exit(-1)
#     if args.protid_file: 
         
#         for one_protid in args.protid:
#             try:
#                 PDBmapper(one_protid, interfacesDB_subset, VCF_subset,
#                         args.out)
#             except IOError:
#                 print('ERROR: Ensembl protein id provided is no supported.')
#                 exit(-1)

##########################
# execute main function  #
##########################
if __name__ == '__main__':
    # measure execution time
    start = timer()
    main()
    end = timer()
    finish = end - start
    print('Congratulations!. PDBmapper has run in', finish, 'seconds.')
